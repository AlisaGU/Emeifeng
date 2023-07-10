#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
gffPath=""
genomePath=""
seq_for_each_cds=""
cdsSeq_for_each_mRNA=""
pepSeq_for_each_mRNA=""
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
seqkit subseq --gtf $gffPath --feature cds $genomePath >$seq_for_each_cds
s1.s1_get_longest_transcript_seq.R $gffPath $seq_for_each_cds $cdsSeq_for_each_mRNA $pepSeq_for_each_mRNA