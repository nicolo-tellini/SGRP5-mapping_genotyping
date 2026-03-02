## 2) function to merge all the gvcf according to stats recorded.
## run it manually ##

source ./config

cd $gvcfstatdir
# create files
> $gvcfstatdir/all.geno.txt
> $gvcfstatdir/df_keep.txt
> $gvcfstatdir/df_discard.txt
# remove files
rm $gvcfstatdir/all.geno.txt $gvcfstatdir/df_keep.txt $gvcfstatdir/df_discard.txt
# add all samples stat on single file 
cat * > $gvcfstatdir/all.geno.txt
# assess strains based on cutoff decided in  the config file, note this are indicative and mark samples as keep and discard (this step require manual curation) 
/usr/bin/time -v Rscript "$basedir/scr/parserstats.r" "$gvcfstatdir" > "$basedir/scr/log/parserstats.log" 2> "$basedir/scr/log/parserstats.err"
# create n files containing each 100 sample names to be merged 
cut -f1 $gvcfstatdir/df_keep.txt | sed "s+$+.$ref.bcftools.filt.plink.gvcf.gz+" | split -l 100 - $gvcfdir/batch_
# enter the gVCF dir
cd $gvcfdir
# define the function for merging gVCF
merge_gvcf() { batch=$1;   bcftools merge -l "$batch" -W -Oz1 -o "$batch.gvcf.gz";}
# export the function
export -f merge_gvcf
# apply the function in a parallel manner so thatt each batch of strain will be merged in a parallel fashion 
ls $gvcfdir/batch_* | parallel -j 4 merge_gvcf

wait
# merge etch batch.gvcf.gz in a single gVCF multi sample conteining all samples 
bcftools merge $gvcfdir/batch_*gz -W -Oz1 -o $gvcfdir/allbatches.gvcf.gz --force-single

wait
# rem interemediate files 
rm $gvcfdir/batch_*
# keep only biallelic variant sites with a frequence of missingness to 0
bcftools view -m2 -M2 -v snps -i 'F_MISSING=0' $gvcfdir/allbatches.gvcf.gz -Oz1 -o $gvcfdir/allbatches.vcf.gz -W
# covert VCF to fasta file 
vcf2phy -i $gvcfdir/allbatches.vcf.gz -p -f -m 10000 -o $gvcfdir/allbatches.fasta
# create a nj tree, this is used for rapid assessment of groups
FastTree -nt -nj $gvcfdir/*.fasta  > $basedir/tree/allbathches.fasttree.nwk 2> "$basedir/scr/log/fasttree.err"
