#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

def parse_arguments():
  """
  Parse command line arguments for variant calling and annotation pipeline.
  
  Returns:
      argparse.Namespace: Parsed command line arguments including:
          - strand: Which DNA strand to consider
          - min_alt_count: Minimum alternative base count threshold
          - min_count: Minimum read depth threshold
          - min_freq: Minimum variant frequency threshold
          - ref: Reference FASTA file path
          - annot: GFF annotation file path
          - prot_attr: Protein name attribute in GFF
          - labels: Sample labels
          - out: Output file path
          - base_files: Input base count files
          - indel_files: Input indel files
  """
  parser = argparse.ArgumentParser(description='Process CSV files to generate variant calls and annotations.')
  parser.add_argument('--strand', choices=['both', 'forward', 'reverse'], default='both', help='Strand to consider: both, forward, or reverse')
  parser.add_argument('--min_alt_count', type=int, default=1, help='Minimum number of counts to consider an alternative base for calling')
  parser.add_argument('--min_count', type=int, default=1, help='Minimum number of reads to consider a position for calling')
  parser.add_argument('--min_freq', type=float, default=0.0, help='Minimum frequency to call an alternative base')
  parser.add_argument('--min_freq_indel', type=float, default=-1.0, help='Minimum frequency to call an indel (default =min_freq)')
  parser.add_argument('--ref', required=True, help='Reference FASTA file')
  parser.add_argument('--annot', required=True, help='GFF annotation file')
  parser.add_argument('--prot_attr', default='ID', help='Attribute with protein name in the GFF')
  parser.add_argument('--labels', required=True, help='Comma-separated sample labels (e.g., spl1,spl2,spl3)')
  parser.add_argument('--out', required=True, help='Output CSV file')
  parser.add_argument('--out_prot', required=False, help='Output FASTA file for proteins')
  parser.add_argument('--vcf', required=False, help='Output VCF file')
  parser.add_argument('--flank_n', type=int, default = 999, help="Number of flanking base to extract at the end of CDS in case of stop codon loss or fs")
  parser.add_argument('--base_files', nargs='+', help='CSV files for each sample in the same order as labels')
  parser.add_argument('--indel_files', nargs='+', help='Indel CSV files for each sample in the same order as labels')
  return parser.parse_args()

def read_sample_csv(filename, sample_label, strand, min_count, min_alt_count, min_freq):
  """
  Read and process base count CSV files to identify variants.
  
  Args:
      filename (str): Path to input CSV file
      sample_label (str): Sample identifier
      strand (str): Which strand to consider ('both', 'forward', or 'reverse')
      min_count (int): Minimum read depth threshold
      min_alt_count (int): Minimum alternative base count threshold
      min_freq (float): Minimum variant frequency threshold
  
  Returns:
      pandas.DataFrame: Processed variant calls with columns:
          - region: Reference sequence identifier
          - pos: Position in reference
          - ref: Reference base
          - alt: Alternative base
          - D:{sample}: Depth at position
          - F:{sample}: Variant frequency
          - {sample}_called: Boolean indicating if variant passes thresholds
  
  Note:
      Strategy for variant calling:
      1. Loads base counts from CSV
      2. Calculates depth based on specified strand
      3. For each position, considers all possible alternative bases
      4. Applies thresholds to determine called variants
  """
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
    ref_base = row['ref'].upper()
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
  """
  Read and process indel CSV files.
  
  Args:
      filename (str): Path to indel CSV file
      sample_label (str): Sample identifier
      depths (dict): Pre-calculated depths per position
      strand (str): Which strand to consider
      min_count (int): Minimum read depth threshold
      min_alt_count (int): Minimum alternative count threshold
      min_freq (float): Minimum variant frequency threshold
  
  Returns:
      pandas.DataFrame: Processed indel calls with depth and frequency information
  """
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
    df.loc[index, "ref"] = row["ref"].upper()

  df[f_col] = df["alt_depth"] / df[d_col]
  df[c_col] = False
  for index, row in df.iterrows():
    if (row[f_col] > min_freq) and (row[d_col] > min_count) and (row['alt_depth'] > min_alt_count):
      df.loc[index, c_col] = True
  df.drop(columns=["forward","reverse","alt_depth"], inplace = True)
  return df

def merge_samples(sample_dfs, labels):
  """
  Merge variant calls from multiple samples.
  
  Args:
      sample_dfs (list): List of DataFrames containing variant calls
      labels (list): Sample labels
  
  Returns:
      pandas.DataFrame: Merged variant calls across all samples
      
  Note:
      Strategy:
      1. Performs outer merge on region, position, reference and alternative bases
      2. Fills missing values appropriately
      3. Keeps only variants called in at least one sample
  """
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
  """
  Parse GFF attribute string into ID and Parent values.
  
  Args:
      attributes (str): GFF9 format attribute string
  
  Returns:
      tuple: (ID, Parent) extracted from attributes
  """
  attr_dict = {}
  for attr in attributes.split(';'):
    if '=' in attr:
      key, value = attr.split('=', 1)
      value = value.split(" ")[0]
      attr_dict[key.strip()] = value.strip()
  return attr_dict.get('ID', None), attr_dict.get('Parent', None)

def parse_gff_to_dataframe(gff_filename):
  # Read the GFF file into a DataFrame
  # How to ignore any additional columns (more than 9) and ensure the first column is reads and not used as index
  gff_df = pd.read_csv(
      gff_filename,
      sep='\t',
      comment='#',
      usecols=range(9),
      header=None,
      names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'],
      dtype={'seqid': str, 'source': str, 'type': str, 'start': int, 'end': int, 'score': str, 'strand': str, 'phase': str, 'attributes': str},
      index_col=False
  )  

  gff_df[['ID', 'Parent']] = gff_df['attributes'].apply(lambda x: pd.Series(extract_attributes(x)))

  return gff_df

def build_cds_dict(gff_df, prot_attr, ref_sequences, flank_n, prot_filename = None):
  """
  Build a nested dictionary of CDS information from GFF annotations.
  
  Args:
      gff_df (pandas.DataFrame): GFF annotation data
      prot_attr (str): Attribute containing protein names
      ref_sequences (dict): Reference sequences
      flank_n (int): Number of flanking bases to include
      prot_filename (str): Path to write proteins
  
  Returns:
      dict: Nested dictionary with structure:
          {sequence_id: {parent_id: {
              'strand': strand,
              'prot_name': protein_name,
              'parts': [[start, end], ...],
              'seq': concatenated_sequence
          }}}
  
  Note:
      Strategy:
      1. Groups CDS features by parent
      2. Validates strand consistency
      3. Orders CDS parts by position considering strand
      4. Concatenates sequences including flanking regions
      5. Handles both forward and reverse strand features
  """

  if prot_filename:
    with open(prot_filename, 'w') as f:
      pass

  # Filter for CDS features
  cds_df = gff_df[gff_df['type'].str.lower() == 'cds']

  # Initialize nested dictionary for CDS features
  cds_dict = defaultdict(lambda: defaultdict(lambda: {'strand': None, 'prot_name': "",'parts': [], 'seqs': [], 'seq': None}))

  for _, row in cds_df.iterrows():
    seqid = row['seqid']
    parent = row['Parent']
    strand = row['strand']
    start = row['start']
    end = row['end']

    # POTENTIAL ISSUE: Complex hierarchical structures in GFF might lead to missing protein names
    # Traverse up the parent hierarchy until we find the protein name
    # This assumes the protein name will be found in an ancestor feature
    prot_name = None
    current_parent = parent
    prot_name = row['attributes'].split(f'{prot_attr}=')[-1].split(';')[0] if f'{prot_attr}=' in row['attributes'] else None
    while current_parent and not prot_name:
      parent_row = gff_df[gff_df['ID'] == current_parent]
      if not parent_row.empty:
        # Extract protein name from attributes
        # POTENTIAL ISSUE: Attribute parsing might fail if format is inconsistent
        prot_name = parent_row.iloc[0]['attributes'].split(f'{prot_attr}=')[-1].split(';')[0] if f'{prot_attr}=' in parent_row.iloc[0]['attributes'] else None
        current_parent = parent_row.iloc[0]['Parent']
      else:
        current_parent = None

    # Store protein name if found
    if prot_name:
      cds_dict[seqid][parent]['prot_name'] = prot_name.split(" ")[0]

    # Validate strand consistency
    # POTENTIAL ISSUE: Multi-exon genes on different strands would raise error
    if parent is None:
      raise ValueError(f"CDS feature at {seqid}:{start}-{end}{strand} has no parent")
    if cds_dict[seqid][parent]['strand'] is None:
      cds_dict[seqid][parent]['strand'] = strand
    elif cds_dict[seqid][parent]['strand'] != strand:
      raise ValueError(f"Inconsistent strand information at {seqid}:{start}-{end}{strand}")

    # Extract and store CDS sequence
    # POTENTIAL ISSUE: Reference sequence might not match GFF coordinates
    cds_dict[seqid][parent]['parts'].append([start, end])
    cds_seq = ref_sequences[seqid].seq[(start - 1):end].upper()
    if strand == "-":
      cds_seq = cds_seq.reverse_complement()
    cds_dict[seqid][parent]['seqs'].append(cds_seq)

  # Process each CDS to create complete sequences
  # add Flanking sequence that extend beyond reference bounds to handle framshift
  # or stop codon loss
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
        flank_seq = ref_sequences[seqid].seq[(flank_start - 1):flank_end].upper()
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
        flank_seq = ref_sequences[seqid].seq[(flank_start - 1):flank_end].reverse_complement().upper()
      if prot_filename:
        prot_seq = None
        for seq in cds_dict[seqid][parent]['seqs']:
          if prot_seq is None:
            prot_seq = seq
          else:
            prot_seq += seq
        print(cds_dict[seqid][parent]['prot_name'])
        prot_seq = prot_seq.translate()
        with open(prot_filename, 'a') as f:
          fasta_id = cds_dict[seqid][parent]['prot_name']
          f.write(f'>{fasta_id}\n')
          for i in range(0, len(prot_seq), 60):
            f.write(str(prot_seq[i:i+60]) + '\n')
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
  """
  Annotate a single variant with protein-level information.
  
  Args:
      row (pandas.Series): Variant information
      cds_info (dict): CDS information including sequence and coordinates
      gff_id (str): GFF feature identifier
  
  Returns:
      dict: Annotation information including:
          - Protein name
          - Amino acid position
          - Reference/alternative codons
          - Reference/alternative amino acids
          - Mutation type
          - GFF ID
  
  Note:
      Strategy:
      1. Determines CDS position and codon context
      2. Handles special cases (indels, frameshifts)
      3. Translates reference and alternative sequences
      4. Classifies mutation type
      5. Special handling for stop codon loss
  """
  # Calculate CDS position by iterating through exons
  # POTENTIAL ISSUE: Position calculation might be off for complex splice patterns
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
    raise ValueError("pos not found in CDS")
  
  # Initialize variant type flags
  is_ins = False
  is_del = False
  is_fs = False
  indel_size = 0

  # Calculate codon position (0-based)
  codon_start = ((cds_pos -1) // 3) * 3 # 0-based
  codon_pos = ((cds_pos -1) % 3) # 0-based

  # Extract reference codon sequence
  ref_codon_seq = []
  for i in range(3):
    ref_codon_seq.append(cds_info['seq'][codon_start + i])
  
  # Create global value to handle frameshift in the global context
  pos2add = None
  
  # Handle different variant types (SNV, insertion, deletion)
  alt_codon_seq = [l.lower() for l in ref_codon_seq]
  indel_size = len(row['alt']) - 1
  if row['alt'].startswith('+'): # Insertion
    is_ins = True
    alt_codon_seq.insert(codon_pos + 1, row['alt'][1:])
    pos2add = codon_start + 3
  elif row['alt'].startswith('-'): # Deletion
    is_del = True
    # define the last deletion codon position (start of codon and codon context)
    del_end_cds_pos = cds_pos + indel_size
    del_end_codon_start = ((del_end_cds_pos -1) // 3) * 3 # 0-based
    del_end_codon_pos = ((del_end_cds_pos -1) % 3) # 0-based
    # ref_codon_seq may include more than 1 codon if span across codons 
    ref_codon_seq = []
    for i in range(codon_start, del_end_codon_start + 1, 3):
      for j in range(3):
        ref_codon_seq.append(cds_info['seq'][i + j])
    # alt_codon shall not include deleted base
    alt_codon_seq = []
    for i in range(codon_start, del_end_codon_start + 1, 3):
      for j in range(3):
        pos2add = i + j
        if (pos2add > codon_start + codon_pos) and (pos2add <= del_end_codon_start + del_end_codon_pos):
          alt_codon_seq.append("-")
        else:
          alt_codon_seq.append(cds_info['seq'][pos2add])
  else: # SNV
    alt_codon_seq[codon_pos] = row['alt']

  # Check for frameshift
  if indel_size % 3 != 0:
    is_fs = True

  # Construct the alternative codon
  if not (is_del or is_ins):
    alt_codon_seq[codon_pos] = row['alt']

  # switch from list to str for codon seqs
  alt_codon_seq = ''.join(alt_codon_seq)
  ref_codon_seq = ''.join(ref_codon_seq)

  # Get amino acids
  ref_aa = str(Seq(ref_codon_seq).translate())
  aa_pos = int(1 + (codon_start / 3)) # 1-based     

  # Handle frameshift consequences
  if is_fs:
    fs_alt_codon_seq = alt_codon_seq
    alt_aa = "?"
    while not alt_aa.endswith("*"):
      if pos2add + 3 > len(cds_info['seq']):
        break
      pos2add += 1
      fs_alt_codon_seq += str(cds_info['seq'][pos2add])
      while len(fs_alt_codon_seq.replace('-','')) % 3 != 0:
        pos2add += 1
        fs_alt_codon_seq += str(cds_info['seq'][pos2add])
      alt_aa = str(Seq(fs_alt_codon_seq.replace('-','')).translate())
  else:
    alt_aa = str(Seq(alt_codon_seq.replace('-','')).translate())

  # Classify mutation type
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

  # Handle stop loss specially
  if mut_type == "stop_lost":
    next_codon_start = codon_start + 3
    while next_codon_start + 2 < len(cds_info['seq']) and (not alt_aa.endswith("*")):
      next_codon = []
      for i in range(3):
        next_codon.append(cds_info['seq'][next_codon_start + i])
      alt_aa += str(Seq("".join(next_codon)).translate())
      next_codon_start += 3

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

def annotate_variants(merged_df, cds_dict):
  """
  Annotate all variants with protein-level information.
  
  Args:
      merged_df (pandas.DataFrame): Merged variant calls
      ref_sequences (dict): Reference sequences
      cds_dict (dict): CDS information dictionary
  
  Returns:
      pandas.DataFrame: Annotated variants with protein-level information
      
  Note:
      Processes each variant to determine if it falls within a CDS
      and adds appropriate annotation or marks as noncoding
  """
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

def write_vcf(annotated_df, output_file, reference_id, labels, source="sav_call"):
  """
  Write variants to VCF format.
  
  Args:
      annotated_df (pandas.DataFrame): Annotated variants DataFrame
      output_file (str): Path to output VCF file
      reference_id (str): Reference sequence identifier
      labels (list): Sample labels
      source (str): Source program name
  
  Note:
      Handles special cases:
      - Multiple annotations per variant (overlapping features)
      - Indels (+ and - prefixes in alt field)
      - Missing values
  """
  from datetime import datetime
  
  # Open output file
  with open(output_file, 'w') as f:
    # Write header
    f.write("##fileformat=VCFv4.2\n")
    f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
    f.write(f"##source={source}\n")
    f.write(f"##reference={reference_id}\n")
    
    # Write INFO and FORMAT definitions
    f.write(
      '##INFO=<ID=ANN,Number=A,Type=String,Description="Annotation per ALT allele: '
      "'Prot_Name|AA_pos|AA_ref|AA_alt|CODON_ref|CODON_alt|GFF_ID' "
      'Multiple features within same ALT separated by semicolon">\n')
    f.write('##FORMAT=<ID=FQ,Number=A,Type=Float,Description="Frequency of each alternative allele">\n')
    f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth at position (ref + all alt)">\n')
    
    # Write column headers
    columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + labels
    f.write('\t'.join(columns) + '\n')
    
    # Group variants by position to handle multiple annotations
    grouped = annotated_df.groupby(['region', 'pos', 'ref', 'alt'])
    
    # Process each variant
    for (region, pos, ref, alt), group in grouped:
      # Skip non-coding variants if they don't have protein annotation
      if group['mut_type'].iloc[0] == 'noncoding' and group['prot_name'].iloc[0] == '':
        continue
          
      # Handle indels
      vcf_alt = alt
      if alt.startswith('+'):  # insertion
        vcf_alt = ref + alt[1:]
      elif alt.startswith('-'):  # deletion
        vcf_ref = ref + alt[1:]
        vcf_alt = ref
        ref = vcf_ref
        alt = vcf_alt
      
      # Build annotation string
      ann_parts = []
      for _, row in group.iterrows():
        if row['prot_name']:  # Only include if there's protein annotation
          ann = f"{row['prot_name']}|{row['aa_pos']}|{row['ref_aa']}|{row['alt_aa']}|{row['ref_codon']}|{row['alt_codon']}|{row['gff_id']}"
          ann_parts.append(ann)
      
      info = f"ANN={';'.join(ann_parts)}" if ann_parts else "."
      
      # Format genotype fields
      format_str = "FQ:DP"
      sample_fields = []
      for label in labels:
        freq = group[f'F:{label}'].iloc[0]
        depth = group[f'D:{label}'].iloc[0]
        if pd.isna(freq):
          freq = 0.0
        sample_fields.append(f"{freq:.4f}:{int(depth)}")
      
      # Write variant line
      fields = [
          region,                  # CHROM
        str(pos),                  # POS
        ".",                       # ID
        ref,                       # REF
        vcf_alt,                   # ALT
        ".",                       # QUAL
        "PASS",                    # FILTER
        info,                      # INFO
        format_str                 # FORMAT
      ] + sample_fields            # Sample columns
      
      f.write('\t'.join(fields) + '\n')

def main():
  """
  Main function orchestrating the variant calling and annotation pipeline.
  
  Workflow:
  1. Parse command line arguments
  2. Load reference sequences and annotations
  3. Process base count and indel files for each sample
  4. Merge variants across samples
  5. Annotate variants with protein-level information
  6. Output results to CSV
  """
  args = parse_arguments()

  labels = args.labels.split(',')
  if len(labels) != len(args.base_files) or len(labels) != len(args.indel_files):
    sys.exit('The number of labels does not match the number of input files.')

  # Read the reference sequence
  ref_sequences = SeqIO.to_dict(SeqIO.parse(args.ref, 'fasta'))

  # Read GFF annotation file
  gff_df = parse_gff_to_dataframe(args.annot)
  cds_dict = build_cds_dict(gff_df, args.prot_attr, ref_sequences, args.flank_n, args.out_prot)

  min_freq_indel = args.min_freq
  if (args.min_freq_indel > 0):
    min_freq_indel = args.min_freq_indel

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
    df_indel = read_indel_csv(indel_file, label, all_depths[label], args.strand, args.min_count, args.min_alt_count, min_freq_indel)
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

  annotated_df = annotate_variants(sorted_df, cds_dict)
  cols = ["region","pos","ref","alt","prot_name","aa_pos","ref_codon","alt_codon","ref_aa","alt_aa","mut_type"]
  for label in labels:
    cols.append(f'F:{label}')
  for label in labels:
    cols.append(f'D:{label}')
  cols.append('gff_id')

  annotated_df.to_csv(args.out, columns=cols, index=False, sep=';', float_format='%.4f')
  if args.vcf:
    write_vcf(annotated_df, args.vcf, args.ref, labels)

if __name__ == '__main__':
  main()