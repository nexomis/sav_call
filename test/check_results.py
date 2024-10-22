import pandas as pd
import argparse
import sys

# Here is the idea, could you analyse it and make sure the position are good, tell in comment what is 0-based and what is 1-based

def get_alt_seq_from_indel(alt, pos, base_seq):
  # parse alt
  sign = alt[0]
  alt_seq = alt[1:]
  if sign == "-":
    del_size = len(alt_seq)
    return(base_seq[0:pos+1]+base_seq[1+pos+del_size:])
  else:
    return(base_seq[0:pos+1]+alt_seq+base_seq[1+pos:])

def compare_indels(expected_region, expected_pos, expected_alt, found_mutations, reference_sequence, tolerance=10, flank=3):
  """
  Compare indels based on expected and found positions and sequences.
  """
  for _, found in found_mutations.iterrows():
      found_region = found['region']
      found_pos = found['pos']
      found_alt = found['alt']

      # Check if positions are within the tolerance range
      if (abs(expected_pos - found_pos) <= tolerance) and (expected_region == found_region):
          # Determine the range for sequence comparison
          start = min(expected_pos, found_pos) - flank -1
          end = max(expected_pos, found_pos) + flank

          # Extract sequences for comparison
          base_seq = reference_sequence[start:end]

          # parse
          found_alt_seq = get_alt_seq_from_indel(found_alt, found_pos - start - 1, base_seq)
          expected_alt_seq = get_alt_seq_from_indel(expected_alt, expected_pos - start - 1, base_seq)

          # Check if the sequences match
          if expected_alt_seq == found_alt_seq:
              return found

  return None


def compare_mutations(expected_indel_file, expected_snp_file, found_file, reference_file):
  # Load expected mutations
  expected_indels = pd.read_csv(expected_indel_file)
  expected_snps = pd.read_csv(expected_snp_file)

  # Load found mutations
  found_mutations = pd.read_csv(found_file, sep=";")

  # Load reference sequence
  reference_sequence = ""
  with open(reference_file, 'r') as ref_file:
    for line in ref_file.readlines():
      line = line.rstrip()
      if line[0] != ">": reference_sequence += line

  # Adjust positions for comparison (0-based to 1-based)
  expected_snps['position'] += 1

  # Prepare a DataFrame for results
  results = []

  # Compare SNPs
  for _, expected in expected_snps.iterrows():
      region = expected['region']
      pos = expected['position']
      alt_base = expected['alt_base']
      expected_prop = expected['proportion']

      # Find matching found mutation
      found = found_mutations[(found_mutations['region'] == region) & (found_mutations['pos'] == pos) & (found_mutations['alt'] == alt_base)]
      if not found.empty:
          found_freq = found['freq'].values[0]
          results.append({
              'region': found['region'].values[0],
              'found_pos': found['pos'].values[0],
              'expected_pos': pos,
              'expected_alt': alt_base,
              'expected_prop': expected_prop,
              'found_alt': found['alt'].values[0],
              'found_freq': found_freq
          })
      else:
          results.append({
              'region': region,
              'found_pos': None,
              'expected_pos': pos,
              'expected_alt': alt_base,
              'expected_prop': expected_prop,
              'found_alt': None,
              'found_freq': 0
          })

  # Compare Indels
  for _, expected in expected_indels.iterrows():
      region = expected['region']
      pos = expected['position']
      alt_base = expected['sequence']
      if expected['type'] == "insertion": 
        alt_base = "+" + expected['sequence']
      if expected['type'] == "deletion": 
        alt_base = "-" + expected['sequence'] 
      
      expected_prop = expected['proportion']
      expected_pos = expected['position']

      # Find matching found mutation
      found = compare_indels(region, pos, alt_base, found_mutations, reference_sequence)
      if found is not None:
          found_freq = found['freq']
          results.append({
              'region': found['region'],
              'found_pos': found['pos'],
              'expected_pos': expected_pos,
              'expected_alt': alt_base,
              'expected_prop': expected_prop,
              'found_alt': found['alt'],
              'found_freq': found_freq
          })
      else:
          results.append({
              'region': region,
              'found_pos': pos,
              'expected_pos': expected_pos,
              'expected_alt': alt_base,
              'expected_prop': expected_prop,
              'found_alt': None,
              'found_freq': 0
          })

  # Track unmatched found mutations
  matched_positions = set((row['found_pos'], row['found_alt']) for row in results if row['found_alt'] is not None)
  unmatched_found = found_mutations[~found_mutations.apply(lambda row: (row['pos'], row['alt']) in matched_positions, axis=1)]

  # Add unmatched found mutations to results
  for _, found in unmatched_found.iterrows():
      results.append({
          'region': found['region'],
          'found_pos': found['pos'],
          'expected_pos': None,
          'expected_alt': None,
          'expected_prop': 0,
          'found_alt': found['alt'],
          'found_freq': found['freq']
      })
    
  # Create a DataFrame from results
  results_df = pd.DataFrame(results)
  results_df['diff_prop'] = results_df['expected_prop'] - results_df['found_freq']
  results_df.sort_values(by='found_pos', ascending=True, inplace=True)

  # Convert expected_prop and found_freq to percentage with 2 decimal places
  results_df['expected_prop'] = results_df['expected_prop'].apply(lambda x: round(x * 100, 2) if pd.notnull(x) else x)
  results_df['found_freq'] = results_df['found_freq'].apply(lambda x: round(x * 100, 2) if pd.notnull(x) else x)
  results_df['diff_prop'] = results_df['diff_prop'].apply(lambda x: round(x * 100, 2) if pd.notnull(x) else x)
  # Output to stdout and stderr based on diff_prop
  stdout_df = results_df[results_df['diff_prop'].abs() < 2]
  stderr_df = results_df[results_df['diff_prop'].abs() >= 2]

  # Output to stdout
  stdout_df.to_csv(sys.stdout, sep='\t', index=False)

  # Output to stderr
  stderr_df.to_csv(sys.stderr, sep='\t', index=False)


def main():
  parser = argparse.ArgumentParser(description='Compare expected mutations with found mutations.')
  parser.add_argument('expected_indel_file', type=str, help='Path to the expected indel mutations file (CSV)')
  parser.add_argument('expected_snp_file', type=str, help='Path to the expected SNP mutations file (CSV)')
  parser.add_argument('found_file', type=str, help='Path to the found mutations file (CSV)')
  parser.add_argument('reference_file', type=str, help='Path to the reference sequence file (FASTA)')

  args = parser.parse_args()

  compare_mutations(args.expected_indel_file, args.expected_snp_file, args.found_file, args.reference_file)

if __name__ == '__main__':
  main()