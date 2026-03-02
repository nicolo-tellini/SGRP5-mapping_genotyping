# SGRP5-genomics
Genomics pipelines
## mapping and calling 
mamba create -n mapping_calling \
    bwa=0.7.17 \
    samtools=1.20 \
    gawk=5.3.1 \
    parallel=20240722 \
    -c conda-forge -c bioconda
