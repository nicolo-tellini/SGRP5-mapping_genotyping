### Set in variable so that it corrrsponds to your gVCF ###
####### MUST BE PRESENT :
# med.gvcf.gz.csi : indice gvcf
# med.gvcf.gz: gVCF
# core.bed : BED of core region
# sample-origin.txt : sample tab country tab lat tab lon
# scripts : dir with the scripts 
in=med.gvcf.gz

cd scripts
for i in *; do
    ln -s scripts/$i ../$i
done
cd ..

###################################################
######      QUAL, FILT, LD AND PRUNING       ######
###################################################

## Extract biallelic SNPs in core genomic regions.
## - COUNT(GT~"1")>0: retain sites where ALT allele 1 is observed in at least one genotype.
## - -m2 -M2: keep  biallelic variants.
## - -R core.bed: restrict analysis to core regions (no subtelomeric regions).
bcftools view -i 'COUNT(GT~"1")>0' -m2 -M2 -R core.bed $in -Oz -o varonly.core.vcf.gz

## Get inbreeding coefficient
plink2 --vcf varonly.core.vcf.gz --het --out het

## Get rate of missing sample   
plink2 --vcf varonly.core.vcf.gz --missing --out varonly.core.missing
rm varonly.core.missing.vmiss
Rscript sample-missing.r

# Filter on missingness.
# --mind 0.1: remove strains with >2.5% missing sites.
plink2 --vcf varonly.core.vcf.gz --mind 0.025 --make-bed --out varonly.core.mind0025

## Get rate of missing var   
plink2 --bfile varonly.core.mind0025 --missing --out varonly.core.mind0025
rm varonly.core.mind0025.smiss
Rscript var-missing.r

# Filter on missingness.
# --geno 0.05: remove variants with >5% missing sites.
plink2 --vcf varonly.core.vcf.gz --geno 0.05 --recode vcf --out varonly.core.mind0025geno005

# Take MAC: add AC field and subset according to the value
bcftools +fill-tags varonly.core.mind0025geno005.vcf -- -t AC > varonly.core.mind0025geno005.mac.vcf
bcftools view -i 'INFO/AC=1' -H varonly.core.mind0025geno005.mac.vcf | wc -l > mac.txt
bcftools view -i 'INFO/AC>1' -H varonly.core.mind0025geno005.mac.vcf | wc -l >> mac.txt

# Compute minor allele frequency (MAF) and take sample numbers 
bcftools query -l varonly.core.mind0025geno005.vcf > samples.txt
bcftools +fill-tags varonly.core.mind0025geno005.vcf -- -t MAF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%MAF\n' > maf.txt
Rscript maf_plor.r

# Convert VCF to PLINK binary format.
# --set-all-var-ids @:# assigns variant IDs as CHR:POS.
plink2 --vcf varonly.core.mind0025geno005.vcf --allow-extra-chr --make-bed --out step0 --set-all-var-ids @:#

# Retain variants with MAC >=2
# MAC 2 for removing singletons
plink2 --bfile step0 --mac 2 --make-bed --out step2

# Compute pairwise LD (r²) within 10 kb windows.
# --ld-window 99999 removes SNP-count limit within the kb window.
plink --bfile step2 --r2 --ld-window 99999 --ld-window-kb 100 --ld-window-r2 0 --out ld

## plot LD
Rscript ld_plot.r

## convert back to VCF 
plink --bfile step2 --recode vcf --out step2

## Use appropriate values for filtering
## 
plink --bfile step2 --indep-pairwise 20 kb 10 0.1 --out pruned

# Create LD-pruned dataset.
plink2 --bfile step2 --extract pruned.prune.in --make-bed --out step_pca

# Principal Component Analysis (first 10 PCs).
plink2 --bfile step_pca --pca 10 --out pca
Rscript pca.r 

# Export non–LD-pruned SNP set (after missingness filtering only).
# Used for window-based PCA to retain local LD structure and rare/private variation.
plink2 --bfile step2 --export vcf bgz --out step2vcf

pdfunite sample-missing.pdf var-missing.pdf maf_plot.pdf ld_decay.pdf pca.pdf merged.pdf

#########################################
######            WINPCA           ######
#########################################

# Window-based PCA (winpca) on non–LD-pruned dataset.
conda activate winpca
mkdir -p winpca
cp step2vcf.vcf.gz winpca
cd winpca
gzip -d step2vcf.vcf.gz
cp ../../rep/chrs.txt .

# chrs.txt format: CHR START END LABEL
# Sliding window PCA: 20 kb window, 500 bp increment.
while read -r line
do
  chr=$(echo $line | cut -d" " -f1)
  start=$(echo $line | cut -d" " -f2)
  end=$(echo $line | cut -d" " -f3)
  label=$(echo $line | cut -d" " -f4)
  echo $chr:$start-$end $label
  /home/ntellini/tools/winpca/winpca pca step2vcf.vcf \
    $chr:$start-$end $label --window_size 20000 --increment 500
done < chrs.txt
conda deactivate
cp ../winpca_plot.r .
Rscript winpca_plot.r

#########################################
######            PIXY             ######
#########################################
cd ..

# pixy requires invariant sites; regenerate core-restricted dataset without ALT filtering.
bcftools view -R core.bed $in -Oz -o pixy.core.gvcf.gz --threads 20
bcftools index pixy.core.gvcf.gz --threads 10

# Define populations.
# Each strain is treated as its own population (enables all-vs-all dXY).
# For FST, biologically meaningful multi-sample groups must be defined.
bcftools query -l pixy.core.gvcf.gz | awk '{print $1"\tpop1"}' > populations.txt
conda activate pixy
pixy --vcf pixy.core.gvcf.gz \
     --populations populations.txt \
     --stats pi watterson_theta tajima_d \
     --window_size 2000 \
     --n_cores 10 \
     --output_folder pixy_out \
     --output_prefix pixy_2kb
conda deactivate
cp -r  scripts/pixy_plot.r pixy_out
cd pixy_out
Rscript pixy_plot.r
cd ..

#########################################
######            IQTREE           ######
#########################################

# Export LD-pruned SNPs to bgzipped VCF for phylogenetic analyses.
plink2 --bfile step_pca --recode vcf bgz --out iqtree_snps
bcftools index iqtree_snps.vcf.gz

# Convert VCF to PHYLIP format for IQ-TREE.
vcf2phy -i iqtree_snps.vcf.gz --output-prefix iqtree_snps

# Maximum-likelihood phylogeny with ascertainment bias correction for SNP-only data.
# -m MFP+ASC: model selection with Lewis correction.
# -B 1000: ultrafast bootstrap.
# -alrt 1000: SH-aLRT branch test.
iqtree2 -s iqtree_snps.min4.phy -m MFP+ASC -B 1000 -alrt 1000 -T 44 --prefix iqtree

# Alternative alignment (e.g., variable sites only).
iqtree2 -s iqtree.varsites.phy -m MFP+ASC -B 1000 -alrt 1000 -T 44 --prefix iqtree2 > iqtree2.log 2> iqtree2.err &
