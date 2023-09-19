# Emeifeng
just for small but engaging tricks, like Emeifeng.

# Tricks
## 1. get cds and pep sequence for each mRNA
> code directory: get_cds_and_pep_sequence_for_each_mRNA


```
chmod +x c1.main.sh s1.get_transcript_seq.R
./c1.main.sh gffPath genomePath FileName_seq_for_each_cds FileName_cdsSeq_for_each_mRNA FileName_pepSeq_for_each_mRNA
```
## 2. compute sfs for [FitCoal](https://www.science.org/doi/10.1126/science.abq7487)
Code directory: sfs_for_fitcoal

- test.vcf.gz: an example of the input file
- sfs_computation/missing_summary.txt: missing site summary for the input vcf
- sfs_computation/0.pdf: sfs distribution
- sfs_computation/unfolded_sfs.txt: unfolded sfs
  
```
# create a directory, named sfs_computation, to place files.

./compute_sfs_for_fitcoal.sh `pwd`/test.vcf.gz "outgroup_shan" `pwd`
```
