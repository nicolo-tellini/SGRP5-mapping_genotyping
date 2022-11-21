# SGRP5-genomics
Genomics pipelines

<p align="center">
  <img src="https://github.com/nicolo-tellini/SGRP5-genomics/blob/main/SGRP5org.png" alt="SGRP5org"/>
</p>

- fastQC evalution
  - fastQC [v.0.11.9](https://github.com/s-andrews/FastQC/releases/tag/v0.11.9)
  - MultiQC [v.1.12](https://github.com/ewels/MultiQC/releases/tag/v1.12)
- *spp*. compositions [sppComp-SGRP5](https://github.com/nicolo-tellini/sppComb) 
- *spp*. repartition 

- mapp. ref
  - bwa mem [v.0.7.17-r1198-dirty](https://github.com/lh3/bwa/releases/tag/v0.7.17)
  - samtools [v.1.14] (fixmate,sort,markdup,view,index) - [workflow FASTQtoBAM](http://www.htslib.org/workflow/fastq.html)
   #### CNVs
  - coverage stats
  - per-chromosome average data
   #### per-sample gVCF
  - bcftools [v.1.15.1](https://github.com/samtools/bcftools/releases/tag/1.15.1) (mpileup,call,view,annotate,merge) 
  - bgzip/tabix v.1.11 
  - PLINK [v1.90b6.21](https://manpages.ubuntu.com/manpages/jammy/man1/plink1.9.1.html)
   ####  Phylogeny
  - R (v.4.0.2) packages: 
    - SNPRelate [v.1.22.0](https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html)
    - gdsfmt [v.1.24.1](https://www.bioconductor.org/packages/release/bioc/html/gdsfmt.html)
    - ape [v.5.6-2](https://cran.r-project.org/web/packages/ape/index.html)
    - treeio [v.1.12.0](https://bioconductor.org/packages/release/bioc/html/treeio.html#:~:text='treeio'%20is%20an%20R%20package,from%20different%20sources%20to%20phylogeny.)
    - ggtree [v.2.2.4](https://bioconductor.org/packages/release/bioc/html/ggtree.html)
