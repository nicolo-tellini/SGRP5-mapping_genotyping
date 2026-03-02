#!/bin/bash

source ./config

ind=$1

genome=$(find "$(realpath $basedir/rep)" -name "*.fa")
 
ref=$(echo $genome | rev | cut -d"/" -f1 | rev | cut -d"." -f1)

bam_file="$mapdir/$ind.$ref.bam"

cd "$gvcfdir"

gvcf_file="$gvcfdir/$ind.$ref.bcftools.gvcf.gz"

gbcf_file="$gvcfdir/$ind.$ref.bcftools.gbcf.gz"

filt_gvcf_file="$gvcfdir/$ind.$ref.bcftools.filt.gvcf.gz"

## gVCF parsing 
bcftools head "$gvcf_file" | grep "##" > "$gvcfdir/head1_$ind"

bcftools head "$gvcf_file" | grep -v "##" > "$gvcfdir/head2_$ind"
# remove path from the sample name
sed -i "s+$mapdir/++g" "$gvcfdir/head2_$ind"
# remove the suffix bam from sample name
sed -i "s+.$ref.bam++g" "$gvcfdir/head2_$ind"
# concatenate the 2 parts of the header 
cat "$gvcfdir/head1_$ind" "$gvcfdir/head2_$ind" > "$gvcfdir/head3_$ind"
# clean
rm "$gvcfdir/head2_$ind" "$gvcfdir/head1_$ind"
# update the header of the gVCF
 bcftools reheader -h "$gvcfdir/head3_$ind" "$gvcf_file" -o $gbcf_file
# index the gBCF
 bcftools index $gbcf_file
# clean
rm "$gvcfdir/head3_$ind" $gvcf_file $gvcf_file.csi
 
## gVCF filtering
bcftools filter -e "STRLEN(REF) > 1 || STRLEN(ALT) > 1" "$gbcf_file" | bcftools filter -e "FORMAT/DP < ${minDP} || QUAL < ${minQUAL}" --set-GTs . |  bcftools annotate -x ID,INFO,FILTER -W -Oz -o "$filt_gvcf_file"
