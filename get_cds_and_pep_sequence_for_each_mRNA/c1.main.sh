#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
gffPath=$1
genomePath=$2
seq_for_each_cds=$3
cdsSeq_for_each_mRNA=$4
pepSeq_for_each_mRNA=$5
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
seqkit subseq --gtf $gffPath --feature cds $genomePath >$seq_for_each_cds
s1.get_transcript_seq.R $gffPath $seq_for_each_cds $cdsSeq_for_each_mRNA $pepSeq_for_each_mRNA
