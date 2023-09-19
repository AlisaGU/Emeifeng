#!/bin/bash

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
vcf=$1
outgroup_name=$2
code_dir=$3
# VARIABLE NAMING for test module TODO:
# vcf="/picb/evolgen/users/gushanshan/projects/sfs_for_fitcoal/test.vcf"
# outgroup_name="shanyang-1_FDSW202130148-1r"
# PROCESS TODO:
## Warning ##
echo "####"
echo "Caution: variants in the input vcf should be biallelic SNPs and from noncoding region"\!
echo "####"

## prepare intermediate files ##

echo $(date)" >>> Step1: prepare intermediate files."
cd $(dirname $vcf)
mkdir -p sfs_computation
cd sfs_computation

if [[ $vcf =~ \.gz$ ]]; then
    zgrep -v "^#" $vcf | awk 'BEGIN{
    OFS="\t"
}{
    printf("%s%s",$1,OFS);
    printf("%s%s",$2,OFS);
    printf("%s%s",$4,OFS);
    printf("%s%s",$5,OFS);

    for(i=10;i<=NF;i++){split($i,a,":");split(a[1],b,"/");printf("%s%s%s%s",b[1],OFS,b[2],OFS)}
    printf(ORS)
}' | gzip >simplified.vcf.gz
else
    grep -v "^#" $vcf | awk 'BEGIN{
    OFS="\t"
}{
    printf("%s%s",$1,OFS);
    printf("%s%s",$2,OFS);
    printf("%s%s",$4,OFS);
    printf("%s%s",$5,OFS);

    for(i=10;i<=NF;i++){split($i,a,":");split(a[1],b,"/");printf("%s%s%s%s",b[1],OFS,b[2],OFS)}
    printf(ORS)
}' | gzip >simplified.vcf.gz
fi

outgroup_index_in_raw_vcf=$(zgrep -v "^##" $vcf | head -n 1 | tr "\t" "\n" | grep -n $outgroup_name | awk -F":" '{print $1}')
outgroup_index_in_simplified_vcf=$(echo 2"*"$outgroup_index_in_raw_vcf"-"15 | bc)
outgroup_discard_index=$(echo $outgroup_index_in_simplified_vcf"+1" | bc)

zcat simplified.vcf.gz | awk -v outgroup=$outgroup_index_in_simplified_vcf -v outgroup_discard_index=$outgroup_discard_index '$outgroup==$outgroup_discard_index' | awk -v outgroup=$outgroup_index_in_simplified_vcf -v outgroup_discard_index=$outgroup_discard_index 'BEGIN{
    OFS="\t"
}$outgroup_discard_index="";$outgroup==0{
    $outgroup="";
    print $0
} $outgroup==1{
    $outgroup="";
    a=$3;
    $3=$4;
    $4=a;
    for(i=5;i<=NF;i++){
        if($i==0){
            $i=1
        }else if($i==1){
            $i=0
        }
    }
    print $0
}' | tr -s '\t' '\t' | gzip >simplified.ancestor.vcf.gz

## compute SFS
echo $(date)" >>> Step2: compute and plot SFS."
$code_dir/s1_compute_sfs.R $(pwd)/simplified.ancestor.vcf.gz

## clean up
echo $(date)" >>> Step3: clean up intermediate files."
rm -rf simplified.ancestor.vcf.gz simplified.vcf.gz
