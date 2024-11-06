#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

def parse_arguments():
  parser = argparse.ArgumentParser(description='Process CSV files to generate variant calls and annotations.')
  parser.add_argument('--strand', choices=['both', 'forward', 'reverse'], default='both', help='Strand to consider: both, forward, or reverse')
  parser.add_argument('--min_alt_count', type=int, default=1, help='Minimum number of counts to consider an alternative base for calling')
  parser.add_argument('--min_count', type=int, default=1, help='Minimum number of reads to consider a position for calling')
  parser.add_argument('--min_freq', type=float, default=0.0, help='Minimum frequency to call an alternative base')
  parser.add_argument('--ref', required=True, help='Reference FASTA file')
  parser.add_argument('--annot', required=True, help='GFF annotation file')
  parser.add_argument('--prot_attr', default='ID', help='Attribute with protein name in the GFF')
  parser.add_argument('--labels', required=True, help='Comma-separated sample labels (e.g., spl1,spl2,spl3)')
  parser.add_argument('--out', required=True, help='Output CSV file')
  parser.add_argument('--flank_n', type=int, default = 99, help="Number of flanking base to extract at the end of CDS in case of stop codon loss or fs")
  parser.add_argument('input_files', nargs='+', help='CSV files for each sample in the same order as labels')
  return parser.parse_args()

def read_sample_csv(filename, sample_label, strand, min_count, min_alt_count, min_freq):
  df = pd.read_csv(filename, sep=';', dtype={'region':str, 'pos':int, 'ref':str})
  
  # Determine depth and base counts based on the strand
  if strand == 'both':
    df['depth_considered'] = df['depth']
    bases = ['A', 'T', 'C', 'G', 'N']
    df_counts = df[['region', 'pos', 'ref', 'depth_considered'] + bases + [b.lower() for b in bases]]
    for base in bases:
      df[base + '_count'] = df[base] + df[base.lower()]
  elif strand == 'forward':
    df['depth_considered'] = df['depth_fw']
    bases = ['A', 'T', 'C', 'G', 'N']
    df_counts = df[['region', 'pos', 'ref', 'depth_considered'] + bases]
    for base in bases:
      df[base + '_count'] = df[base]
  elif strand == 'reverse':
    df['depth_considered'] = df['depth_rv']
    bases = [b.lower() for b in ['A', 'T', 'C', 'G', 'N']]
    df_counts = df[['region', 'pos', 'ref', 'depth_considered'] + bases]
    for base in ['A', 'T', 'C', 'G', 'N']:
      df[base + '_count'] = df[base.lower()]
  else:
      sys.exit('Invalid strand option.')
  
  # Generate alternative alleles

  df_region = []
  df_pos = []
  df_ref = []
  df_alt = []
  df_depth = []
  df_freq = []
  df_called = []
  for index, row in df.iterrows():
    total_depth = row['depth_considered']
    ref_base = row['ref']
    for alt_base in ['A', 'T', 'C', 'G', 'N']:
      if alt_base != ref_base:
        alt_count = row[alt_base + '_count']
        freq = alt_count / total_depth if total_depth > 0 else None
        if (freq is None):
          called = False
        else:
          called = (freq >= min_freq) and (alt_count >= min_alt_count) and (total_depth > min_count)
        df_region.append(row['region'])
        df_pos.append(row['pos'])
        df_ref.append(ref_base)
        df_alt.append(alt_base)
        df_depth.append(total_depth)
        df_freq.append(freq)
        df_called.append(called)
        
  return(pd.DataFrame({
      'region': df_region,
      'pos': df_pos,
      'ref': df_ref,
      'alt': df_alt,
      f'D:{sample_label}': df_depth,
      f'F:{sample_label}': df_freq,
      f'{sample_label}_called': df_called
  }))

def merge_samples(sample_dfs, labels):
  # Merge all sample DataFrames on region, pos, ref, alt
  from functools import reduce

  merged_df = reduce(lambda left, right: pd.merge(left, right, on=['region', 'pos', 'ref', 'alt'], how='outer'), sample_dfs)
  
  # Fill missing values
  cols2drop = ["called_any"]
  for label in labels:
    merged_df[f'D:{label}'] = merged_df.get(f'D:{label}', 0)
    merged_df[f'F:{label}'] = merged_df.get(f'F:{label}', None)
    merged_df[f'{label}_called'] = merged_df.get(f'{label}_called', False)
    cols2drop.append(f'{label}_called')
  
  called_columns = [f'{label}_called' for label in labels]
  merged_df['called_any'] = merged_df[called_columns].any(axis=1)
  merged_df = merged_df[merged_df['called_any']]
  merged_df = merged_df.drop(columns=cols2drop)
  return merged_df

def extract_attributes(attributes):
      attr_dict = {}
      for attr in attributes.split(';'):
          if '=' in attr:
              key, value = attr.split('=', 1)
              attr_dict[key.strip()] = value.strip()
      return attr_dict.get('ID', None), attr_dict.get('Parent', None)

def parse_gff_to_dataframe(gff_filename):
  # Read the GFF file into a DataFrame
  gff_df = pd.read_csv(
      gff_filename,
      sep='\t',
      comment='#',
      header=None,
      names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'],
      dtype={'seqid': str, 'source': str, 'type': str, 'start': int, 'end': int, 'score': str, 'strand': str, 'phase': str, 'attributes': str}
  )  

  gff_df[['ID', 'Parent']] = gff_df['attributes'].apply(lambda x: pd.Series(extract_attributes(x)))

  return gff_df

def build_cds_dict(gff_df, prot_attr, ref_sequences, flank_n):
  # Filter for CDS features
  cds_df = gff_df[gff_df['type'].str.lower() == 'cds']

  # Group CDS by their Parent
  cds_dict = defaultdict(lambda: defaultdict(lambda: {'strand': None, 'prot_name': "",'parts': [], 'seqs': [], 'seq': None}))

  for _, row in cds_df.iterrows():
    seqid = row['seqid']
    parent = row['Parent']
    strand = row['strand']
    start = row['start']
    end = row['end']

    # Find the protein name by looking up the hierarchy
    prot_name = None
    current_parent = parent
    while current_parent and not prot_name:
      parent_row = gff_df[gff_df['ID'] == current_parent]
      if not parent_row.empty:
        prot_name = parent_row.iloc[0]['attributes'].split(f'{prot_attr}=')[-1].split(';')[0] if f'{prot_attr}=' in parent_row.iloc[0]['attributes'] else None
        current_parent = parent_row.iloc[0]['Parent']
      else:
        current_parent = None

    if prot_name:
      cds_dict[seqid][parent]['prot_name'] = prot_name

    # Ensure consistent strand information
    if parent is None:
      raise ValueError(f"CDS feature at {seqid}:{start}-{end}{strand} has no parent")
    if cds_dict[seqid][parent]['strand'] is None:
      cds_dict[seqid][parent]['strand'] = strand
    elif cds_dict[seqid][parent]['strand'] != strand:
      raise ValueError(f"Inconsistent strand information at {seqid}:{start}-{end}{strand}")

    # Add the CDS part
    cds_dict[seqid][parent]['parts'].append([start, end])
    cds_seq = ref_sequences[seqid].seq[(start - 1):end]
    if strand == "-":
      cds_seq = cds_seq.reverse_complement()
    cds_dict[seqid][parent]['seqs'].append(cds_seq)

  # Sort parts based on strand
  flank_seq = None
  for seqid in cds_dict:
    for parent in cds_dict[seqid]:
      if cds_dict[seqid][parent]['strand'] == '+':
        cds_dict[seqid][parent]['parts'], cds_dict[seqid][parent]['seqs'] = zip(
          *sorted(
            zip(cds_dict[seqid][parent]['parts'], cds_dict[seqid][parent]['seqs']),
            key=lambda x: x[0][0]
          )
        )
        flank_start = 1 + max(part[1] for part in cds_dict[seqid][parent]['parts'])
        flank_end = min(flank_start + flank_n - 1, len(ref_sequences[seqid]))
        flank_seq = ref_sequences[seqid].seq[(flank_start - 1):flank_end]
      else:
        cds_dict[seqid][parent]['parts'], cds_dict[seqid][parent]['seqs'] = zip(
          *sorted(
            zip(cds_dict[seqid][parent]['parts'], cds_dict[seqid][parent]['seqs']),
            key=lambda x: x[0][1],
            reverse=True
          )
        )
        flank_end = min(part[0] for part in cds_dict[seqid][parent]['parts']) - 1
        flank_start = max(1, flank_end - flank_n + 1)
        flank_seq = ref_sequences[seqid].seq[(flank_start - 1):flank_end].reverse_complement()
      for seq in cds_dict[seqid][parent]['seqs']:
        if cds_dict[seqid][parent]['seq'] is None:
          cds_dict[seqid][parent]['seq'] = seq
        else:
          cds_dict[seqid][parent]['seq'] += seq
      del cds_dict[seqid][parent]['seqs']
      cds_dict[seqid][parent]['seq'] += flank_seq
      

  # Check for non-standard parent types
  for _, group in cds_df.groupby('Parent'):
    parent_type = gff_df[gff_df['ID'].isin(group['Parent'])]['type'].str.lower().unique()
    if not any(pt in ['mrna', 'transcript'] for pt in parent_type):
      print(f"Warning: CDS with parent type(s) {parent_type} not standard (neither mRNA nor transcript)")

  return cds_dict

def annotate(row, cds_info, gff_id):

  # extract the cds seq by iterating over parts
  cds_offset = 1
  cds_pos = None
  for part in cds_info['parts']:
    if part[0] <= row['pos'] <= part[1]:
      if cds_info['strand'] == "+":
        cds_pos = cds_offset + row['pos'] - part[0]
      else:
        cds_pos = cds_offset + part[1] - row['pos']
    else:
      cds_offset += 1 + part[1] - part[0]

  if cds_pos is None:
    print(cds_info)
    raise ValueError("pos not found in CDS")

  # Determine codon position
  codon_start = ((cds_pos -1) // 3) * 3 # 0-based
  codon_pos = ((cds_pos -1) % 3) # 0-based

  # Get the reference codon
  ref_codon_seq = []
  for i in range(3):
    ref_codon_seq.append(cds_info['seq'][codon_start + i])

  # Construct the alternative codon
  alt_codon_seq = ref_codon_seq.copy()
  alt_codon_seq[codon_pos] = row['alt']
  alt_codon_seq = ''.join(alt_codon_seq)
  ref_codon_seq = ''.join(ref_codon_seq)

  # Get amino acids
  ref_aa = str(Seq(ref_codon_seq).translate())
  alt_aa = str(Seq(alt_codon_seq).translate())
  aa_pos = int(1 + (codon_start / 3)) # 1-based

  # Determine mutation type
  if ref_aa == alt_aa:
    mut_type = 'synonymous'
  elif alt_aa == '*':
    mut_type = 'nonsense'
  elif ref_aa == '*':
    mut_type = 'stop_lost'
  else:
    mut_type = 'nonsynonymous'

  if mut_type == "stop_lost":
    next_codon_start = codon_start + 3
    next_aa = ""
    while next_codon_start + 2 < len(cds_info['seq']) and next_aa != "*":
      next_codon = []
      for i in range(3):
        next_codon.append(cds_info['seq'][next_codon_start + i])
      next_codon = "".join(next_codon)
      next_aa = str(Seq(next_codon).translate())
      alt_aa += next_aa

  # Return the annotated row
  return {
    'prot_name': cds_info['prot_name'],
    'aa_pos': aa_pos,
    'ref_codon': ref_codon_seq,
    'alt_codon': alt_codon_seq,
    'ref_aa': ref_aa,
    'alt_aa': alt_aa,
    'mut_type': mut_type,
    'gff_id': gff_id
  }

def annotate_variants(merged_df, ref_sequences, cds_dict):

  annotated_rows = []

  for _, row in merged_df.iterrows():
    annotated = False
    region = row["region"]
    pos = row["pos"]

    if region in cds_dict:
      for gff_id, cds_info in cds_dict[region].items():
        for part in cds_info['parts']:
          if part[0] <= pos <= part[1]:
            annotated_row = annotate(row, cds_info, gff_id)
            annotated_rows.append({**row, **annotated_row})
            annotated = True
      if not annotated:
        # If not annotated, append the row with empty annotation fields
        annotated_rows.append({**row, 'prot_name': '', 'aa_pos': pd.NA, 'ref_codon': '', 'alt_codon': '', 'ref_aa': '', 'alt_aa': '', 'mut_type': 'noncoding', 'gff_id': ''})

  # Convert the list of annotated rows to a DataFrame
  annotated_df = pd.DataFrame(annotated_rows)
  return annotated_df

def main():
  args = parse_arguments()

  labels = args.labels.split(',')
  if len(labels) != len(args.input_files):
    sys.exit('The number of labels does not match the number of input files.')

  # Read the reference sequence
  ref_sequences = SeqIO.to_dict(SeqIO.parse(args.ref, 'fasta'))

  # Read GFF annotation file
  gff_df = parse_gff_to_dataframe(args.annot)
  cds_dict = build_cds_dict(gff_df, args.prot_attr, ref_sequences, args.flank_n)

  sample_dfs = []
  for filename, label in zip(args.input_files, labels):
    df_sample = read_sample_csv(filename, label, args.strand, args.min_count, args.min_alt_count, args.min_freq)
    # Rename columns to include sample label
    df_sample = df_sample.rename(columns={
      'depth': f'D:{label}',
      'freq': f'F:{label}',
      'called': f'{label}_called'
    })
    sample_dfs.append(df_sample)

  merged_df = merge_samples(sample_dfs, labels)
  annotated_df = annotate_variants(merged_df, ref_sequences, cds_dict)
  cols = ["region","pos","ref","alt","prot_name","aa_pos","ref_codon","alt_codon","ref_aa","alt_aa","mut_type"]
  for label in labels:
    cols.append(f'F:{label}')
  for label in labels:
    cols.append(f'D:{label}')
  cols.append('gff_id')
  annotated_df.to_csv(args.out, columns=cols, index=False, sep=';', float_format='%.4f')

if __name__ == '__main__':
  main()