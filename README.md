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
```
./compute_sfs_for_fitcoal.sh `pwd`/test.vcf.gz "shanyang-1_FDSW202130148-1r" `pwd`
```
