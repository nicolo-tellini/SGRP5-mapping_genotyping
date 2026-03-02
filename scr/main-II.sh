## 2) Merge all gVCFs based on per-sample statistics.
## Execute this section manually.

source ./config

cd $gvcfstatdir

# Initialize empty summary files (will be regenerated below)
> $gvcfstatdir/all.geno.txt
> $gvcfstatdir/df_keep.txt
> $gvcfstatdir/df_discard.txt

# Remove potential residual files from previous runs to avoid duplication
rm $gvcfstatdir/all.geno.txt $gvcfstatdir/df_keep.txt $gvcfstatdir/df_discard.txt

# Remove potential residual files from previous runs to avoid duplication
cat * > $gvcfstatdir/all.geno.txt

# Classify samples according to QC cutoffs defined in config.
# Output files: df_keep.txt and df_discard.txt.
# NOTE: thresholds are indicative and require manual validation before proceeding.
/usr/bin/time -v Rscript "$basedir/scr/parserstats.r" "$gvcfstatdir" > "$basedir/scr/log/parserstats.log" 2> "$basedir/scr/log/parserstats.err"

# Split retained samples into batches of 100 to control memory usage during merging
# Each line contains a full path to the corresponding filtered gVCF
cut -f1 $gvcfstatdir/df_keep.txt | sed "s+$+.$ref.bcftools.filt.plink.gvcf.gz+" | split -l 100 - $gvcfdir/batch_

# Move to gVCF directory for merging
cd $gvcfdir

# Define a function to merge one batch of gVCFs
# -W writes index on the fly
# -Oz1 enables compressed output (level 1 for faster I/O)
merge_gvcf() { batch=$1;   bcftools merge -l "$batch" -W -Oz1 -o "$batch.gvcf.gz";}

# Export the function for GNU parallel
export -f merge_gvcf
# Merge batches in parallel (adjust -j according to available memory and CPUs)
ls $gvcfdir/batch_* | parallel -j 4 merge_gvcf

wait
# Merge all batch-level gVCFs into a final multi-sample gVCF
bcftools merge $gvcfdir/batch_*gz -W -Oz1 -o $gvcfdir/allbatches.gvcf.gz

wait

# Remove intermediate batch files to free disk space
rm $gvcfdir/batch_*

# Retain only biallelic SNPs with no missing data across samples
# This produces a strict dataset suitable for phylogenetic inference
bcftools view -m2 -M2 -v snps -i 'F_MISSING=0' $gvcfdir/allbatches.gvcf.gz -Oz1 -o $gvcfdir/allbatches.vcf.gz -W

# Convert filtered VCF to FASTA alignment
# -m 10000 controls maximum missing sites per sequence before exclusion
python3 $basedir/scr/vcf2phylip.py -i $gvcfdir/allbatches.vcf.gz -p -f -m 10000 -o $gvcfdir/allbatches.fasta

# Build a Neighbor-Joining tree for rapid population structure assessment
# FastTree is used here for speed rather than maximum likelihood optimization
FastTree -nt -nj $gvcfdir/*.fasta  > $basedir/allbathches.fasttree.nwk 2> "$basedir/scr/log/fasttree.err"
