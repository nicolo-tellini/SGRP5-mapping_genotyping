#!/bin/bash

source ./config

ind=$1

genome=$(find "$(realpath $basedir/rep)" -name "*.fa")
 
ref=$(echo $genome | rev | cut -d"/" -f1 | rev | cut -d"." -f1)

bam_file="$mapdir/$ind.$ref.bam"

filt_gvcf_plink_file="$gvcfdir/$ind.$ref.bcftools.filt.plink.gvcf.gz"

files=($filt_gvcf_plink_file $bam_file)

if [[ $(find "${files[@]}" -type f -not -empty 2>/dev/null | wc -l) -eq ${#files[@]} ]]; then
    echo "All files exist. Storing $ind in cps.txt"
    echo "$ind" >> "$cpsdir"/cps.txt
else
    echo "One or more files are missing or empty. NOT storing $ind in cps.txt."
fi
