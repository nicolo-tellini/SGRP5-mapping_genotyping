#!/bin/bash

source ./config

ind=$1

genome=$(find "$(realpath $basedir/rep)" -name "*.fa")
 
ref=$(echo $genome | rev | cut -d"/" -f1 | rev | cut -d"." -f1)

bam_file="$mapdir/$ind.$ref.bam"

cd "$mapdir"

gvcf_file="$gvcfdir/$ind.$ref.bcftools.gvcf.gz"

## Genotype positions BCFTOOLS
bcftools mpileup -E -Ou -q $minMQ -Q $minBQ -a DP -f "$basedir/rep/$ref.genome.fa" "$bam_file" | \
 
bcftools call -mOz1 --ploidy 2 -W -o "$gvcf_file"
