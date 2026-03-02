# SGRP5-genomics

## the env 
```
mamba create -n mapping_calling \
    bwa=0.7.17 \
    samtools=1.20 \
    gawk=5.3.1 \
    parallel=20240722 \
    -c conda-forge -c bioconda
```
## the config file 

```
###########################################################
#                   CONFIGURATION FILE                    #
#                                                         #
# This configuration file sets parameters for sample      #
# processing and analysis. Modify these values according  #
# to your project's requirements.                         #
#                                                         #
# Each parameter is explained briefly to guide            #
# adjustments. Please refer to documentation if further   #
# clarification is needed.				                  #
#                                			              #
# Contact me at nicolo.tellini.2@gmail.com		          #
###########################################################

##########################################################
#		             USER VARIABLE	                     #
##########################################################

# Number of threads for each parallel process		 
nThreads=2

# Number of samples to process in parallel
nSamples=50

# Minimum quality score for variants
minQUAL=20  

# Minimum mapping quality
minMQ=5  

# Ploidy level of the samples
ploidy=2  

# Minimum base quality score
minBQ=30  

# Threshold for heterozygous variants
hetvars_threshold=2500  

# Reference position threshold
ref_pos_threshold=11000000 

## Reference position threshold
## alt_homo_threshold=80000

##########################################################
#    	    DO NOT CHANGE THE VARIABLE BELOW         	 #
##########################################################

basedir=$(dirname "$(pwd)")
mapdir="$basedir/map"
tmpdir="$basedir/tmp"
logdir="$basedir/scr/log"
cpsdir="$basedir/cps"
statdir="$basedir/stat"
gvcfdir="$basedir/gVCF"
seqdir="$basedir/seq"
bamstatdir="$basedir/stat/bam"
gvcfstatdir="$basedir/stat/gvcf"
nasdir="$basedir/nas"
treedir="$basedir/tree"

```
## dir tree
```
.
├── cps
├── gVCF
├── map
├── rep
├── scr
├── seq
├── stat
├── tmp
```
