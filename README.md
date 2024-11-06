# SAV-CALL Strand-aware Variant caller

## Build

```
make
```

## Use

```
Usage: sav_call [options] --bam <input.bam> --reference <reference.fasta>       
                                                                                
  Version 0.2.0                                                                 
                                                                                
Required arguments:                                                             
  --prefix-out <string>                Prefix for output files with called      
                                         variants. 3 files for forward, reverse 
                                         and both strands, respectively:        
                                         - fwd: <prefix>.fwd.csv                
                                         - rev: <prefix>.rev.csv                
                                         - both: <prefix>.both.csv              
  --bam <input.bam>                    Input BAM file                           
  --reference <reference.fasta>        Reference FASTA used for mapping         
                                                                                
Optional arguments:                                                             
  --base-csv <file>                    Base counts and depths for both forward  
                                         and reverse strand with all pileup     
                                         position to CSV file. Position without 
                                         reads are not reported. Gaps in reads  
                                         are accounted in depths                
  --indel-csv <file>                   Output indel counts to CSV file          
  --keep-read-extremities              Keep read extremities (default: false)   
  --keep-duplicate                     Keep duplicate reads (default: false)    
  --keep-secondary                     Keep secondary mapping (default: false)  
  --min-qual <int>                     Minimum base quality to count            
                                         (default: 0)                           
  --R1-strand <forward|reverse>        Set R1 strand (default: forward)         
  --R2-strand <forward|reverse>        Set R2 strand (default: reverse)         
  --bam <input.bam>                    Input BAM file (required)                
  --reference <reference.fasta>        Reference FASTA file (required)          
  --min-freq <float;float;float>       Minimum frequency to call a variant      
                                         (default: 0.02;0.02;0.02)              
                                         fwd;rev;both                           
  --min-count <int;int;int>              Minimum counts to call a variant       
                                         (default: 20;20;20)                    
                                         fwd;rev;both                           
  --min-alt-count <int;int;int>        Minimum alternative counts to call       
                                         (default: 10;10;10)                    
                                         fwd;rev;both                           
  --max-n-pileup <int>                 Maximum reads in pileup                  
                                         (default: 1000000)                     
```

## Example / test

```
make test
```

Coverage report generated in 'test/out' directory.
Check 'test/check*.err' and 'test/check*.err' to see results on sim data.

## Code coverage

```
Overall coverage rate:
  lines......: 79.0% (793 of 1004 lines)
  functions..: 94.5% (397 of 420 functions)
```

## You can use the python script to run a batch analysis with annotation

```
python ../scripts/call.py \
  --labels P0,P1,P2,P3 \
  --strand both \
  --min_alt_count 10 \
  --min_count 20 \
  --min_freq 0.03 \
  --prot_attr gene \
  --ref ref/refseq.fa \
  --annot ref/refseq.gff \
  --out test.python.call.csv \
  --base_files P0.base.csv P1.base.csv P2.base.csv P3.base.csv \
  --indel_files P0.indel.csv P1.indel.csv P2.indel.csv P3.indel.csv
```
