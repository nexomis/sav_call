#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
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
  parser.add_argument('--base_files', nargs='+', help='CSV files for each sample in the same order as labels')
  parser.add_argument('--indel_files', nargs='+', help='Indel CSV files for each sample in the same order as labels')
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
  for _, row in df.iterrows():
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

def read_indel_csv(filename, sample_label, depths, strand, min_count, min_alt_count, min_freq):
  df = pd.read_csv(filename, sep=';', dtype={'region':str, 'pos':int, 'ref':str})
  d_col = f'D:{sample_label}'
  f_col = f'F:{sample_label}'
  c_col = f'{sample_label}_called'
  df[d_col] = 0
  df[f_col] = 0
  if strand == 'both':
    df["alt_depth"] = df["forward"] + df["reverse"]
  elif strand == 'reverse':
    df["alt_depth"] = df["reverse"]
  elif strand == 'forward':
    df["alt_depth"] = df["forward"]
  else: 
    raise ValueError("Wrong strand args")
  #This is not modifying the row in the df, how to do it ?
  for index, row in df.iterrows():
    df.loc[index, d_col] = depths[row["region"]][row["pos"] - 1]

  df[f_col] = df["alt_depth"] / df[d_col]
  df[c_col] = False
  for index, row in df.iterrows():
    if (row[f_col] > min_freq) and (row[d_col] > min_count) and (row['alt_depth'] > min_alt_count):
      df.loc[index, c_col] = True
  df.drop(columns=["forward","reverse","alt_depth"], inplace = True)
  return df

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
  merged_df.drop(columns=cols2drop, inplace = True)
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
  
  is_ins = False
  is_del = False
  is_fs = False
  indel_size = 0

  # Determine codon position
  codon_start = ((cds_pos -1) // 3) * 3 # 0-based
  codon_pos = ((cds_pos -1) % 3) # 0-based

  # Get the reference codon
  ref_codon_seq = []
  for i in range(3):
    ref_codon_seq.append(cds_info['seq'][codon_start + i])
  
  pos2add = None
  
  alt_codon_seq = ref_codon_seq.copy()
  if row['alt'].startswith('+'):
    is_ins = True
    indel_size = len(row['alt']) - 1
    alt_codon_seq.insert(codon_pos, row['alt'][1:])
    pos2add = codon_start + 3
  elif row['alt'].startswith('-'):
    is_del = True
    indel_size = len(row['alt']) - 1
    alt_codon_seq = alt_codon_seq[:codon_pos] + alt_codon_seq[codon_pos + len(row['alt'][1:]):]
  else:
    alt_codon_seq[codon_pos] = row['alt']

  if indel_size % 3 != 0:
    is_fs = True

  # Construct the alternative codon
  if not (is_del or is_ins):
    alt_codon_seq[codon_pos] = row['alt']

  alt_codon_seq = ''.join(alt_codon_seq)
  ref_codon_seq = ''.join(ref_codon_seq)

  # Get amino acids
  ref_aa = str(Seq(ref_codon_seq).translate())
  aa_pos = int(1 + (codon_start / 3)) # 1-based

  # In case of deletion ref_codon shall be built differently
  if is_del:
    del_end_cds_pos = cds_pos + indel_size -1
    del_end_codon_start = ((del_end_cds_pos -1) // 3) * 3 # 0-based
    del_end_codon_pos = ((del_end_cds_pos -1) % 3) # 0-based
    ref_codon_seq = []
    for i in range(codon_start, del_end_codon_start + 1, 3):
      for j in range(3):
        print(cds_info['seq'][i + j])
        ref_codon_seq.append(cds_info['seq'][i + j])
    alt_codon_seq = []
    pos2add = 0
    for j in range(0,codon_pos + 1):
      pos2add = codon_start + j
      alt_codon_seq.append(cds_info['seq'][pos2add])
    if del_end_codon_pos != 2:
      for j in range(del_end_codon_pos + 1,3):
        pos2add = del_end_codon_start + j
        alt_codon_seq.append(cds_info['seq'][pos2add])
    else:
      pos2add = del_end_codon_pos + 2

    alt_codon_seq = ''.join(alt_codon_seq)
    ref_codon_seq = ''.join(ref_codon_seq)
    ref_aa = str(Seq(ref_codon_seq).translate())

  if is_fs:
    fs_alt_codon_seq = alt_codon_seq
    alt_aa = "?"
    while not alt_aa.endswith("*"):
      if pos2add + 3 > len(cds_info['seq']):
        break
      pos2add += 1
      fs_alt_codon_seq += str(cds_info['seq'][pos2add])
      while len(fs_alt_codon_seq) % 3 != 0:
        pos2add += 1
        fs_alt_codon_seq += str(cds_info['seq'][pos2add])
      alt_aa = str(Seq(fs_alt_codon_seq).translate())
  else:
    alt_aa = str(Seq(alt_codon_seq).translate())

  if is_fs:
    mut_type = "frameshift"
  else:
    if is_ins:
      mut_type = "insertion"
    elif is_del:
      mut_type = "deletion" 
    elif ref_aa == alt_aa:
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
  if len(labels) != len(args.base_files) or len(labels) != len(args.indel_files):
    sys.exit('The number of labels does not match the number of input files.')

  # Read the reference sequence
  ref_sequences = SeqIO.to_dict(SeqIO.parse(args.ref, 'fasta'))

  # Read GFF annotation file
  gff_df = parse_gff_to_dataframe(args.annot)
  cds_dict = build_cds_dict(gff_df, args.prot_attr, ref_sequences, args.flank_n)

  sample_dfs = []
  indel_dfs = []
  all_depths = dict()
  for base_file, indel_file, label in zip(args.base_files, args.indel_files ,labels):
    all_depths[label] = dict()
    for region,record in ref_sequences.items():
      all_depths[label][region] = np.zeros(len(record.seq),dtype= np.uint32)
    df_sample = read_sample_csv(base_file, label, args.strand, args.min_count, args.min_alt_count, args.min_freq)
    for _, row in df_sample.iterrows():
      all_depths[label][row["region"]][row["pos"] - 1] = row[f'D:{label}']
    df_indel = read_indel_csv(indel_file, label, all_depths[label], args.strand, args.min_count, args.min_alt_count, args.min_freq)
    sample_dfs.append(df_sample)
    indel_dfs.append(df_indel)

  merged_df = merge_samples(sample_dfs, labels)
  merged_indels = merge_samples(indel_dfs, labels)
  for index, row in merged_indels.iterrows():
    for label in labels:
      d_col = f'D:{label}'
      if pd.isnull(row[d_col]):
        merged_indels.loc[index, d_col] = all_depths[label][row["region"]][row["pos"] - 1]
        merged_indels.loc[index, f'F:{label}'] = 0

  combined_df = pd.concat([merged_df, merged_indels])
  sorted_df = combined_df.sort_values(by=['region', 'pos'], ascending=[True, True])
  # Reset index if needed
  sorted_df = sorted_df.reset_index(drop=True)

  annotated_df = annotate_variants(sorted_df, ref_sequences, cds_dict)
  cols = ["region","pos","ref","alt","prot_name","aa_pos","ref_codon","alt_codon","ref_aa","alt_aa","mut_type"]
  for label in labels:
    cols.append(f'F:{label}')
  for label in labels:
    cols.append(f'D:{label}')
  cols.append('gff_id')
  annotated_df.to_csv(args.out, columns=cols, index=False, sep=';', float_format='%.4f')

if __name__ == '__main__':
  main()