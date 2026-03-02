# SGRP5-genomics

## the env 
```
mamba create -n mapping_calling \
    bwa=0.7.17 \
    samtools=1.20 \
    gawk=5.3.1 \
    bcftools=1.20 \
    parallel=20240722 \
    fasttree=2.2.0 \
    -c conda-forge -c bioconda
```
In addition, https://github.com/edgardomortiz/vcf2phylip is required.
Be sure ```vcf2phylip.py``` is in scr and you can execute it.

## dir tree
```
.
├── rep
├── seq
├── scr
```
