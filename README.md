# SGRP5-genomics

## the env 
```
mamba create -n mapping_calling \
    bwa=0.7.17 \
    samtools=1.20 \
    gawk=5.3.1 \
    bcftools=1.20 \
    parallel=20240722 \
    plink2=2.0.0a.6.9 \
    iqtree=2.4.0 \
    pixy=2.0.0.beta14 \
    fasttree=2.2.0 \
    -c conda-forge -c bioconda
```
In addition, 

- [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip) – tool for VCF to PHYLIP or FASTA
- [winpca](https://github.com/MoritzBlumer/winpca) – window-based PCA along the genome 

## dir tree
```
.
├── rep
├── seq
├── scr
```

Main part I run with bash. 

Main part II and III need to be run manually, interventinons and parameters tuning might be required. 

popGen.sh also requires manul run and flags adjusted according to the population
