# SAV-CALL Stand-aware Variant caller

## Build

```
make
```

## Use

```
Usage: variant_caller [options] --bam <input.bam> --reference <reference.fasta>
Options:
  --count-read-extremities             Count read extremities (default: false)
  --keep-duplicate                     Keep duplicate reads (default: false)
  --keep-secondary                     Keep secondary mapping (default: false)
  --base-csv <file>                    Output base counts to CSV file
  --indel-csv <file>                   Output indel counts to CSV file
  --called-variant-csv <file>          Output called variants to CSV file
  --min-qual <int>                     Minimum base quality to count (default: 0)
  --R1-strand <forward|reverse>        Set R1 strand (default: forward)
  --R2-strand <forward|reverse>        Set R2 strand (default: reverse)
  --bam <input.bam>                    Input BAM file (required)
  --reference <reference.fasta>        Reference FASTA file (required)
  --min-freq <float>                   Minimum frequency to call a variant (default: 0.02)
  --call-strand <forward|reverse|both> Strand to apply thresholds (default: both)
  --min-count <int>                    Minimum count to report in outputs (default: 20)
  --min-alt-count <int>                Minimum alternative counts to call (default: 10)
  --maw-n-pileup <int>                 Maximum reads in pileup (default: 1000000)
```

## Example / test

```
cd test/
bash gene_test_data.sh

```
