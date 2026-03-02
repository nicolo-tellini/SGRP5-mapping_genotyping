#!/bin/bash

set -x 

source ./config

# define path to referrnce genome
genome=$(find "$(realpath $basedir/rep)" -name "*.fa")
# define reference genome name
ref=$(echo $genome | rev | cut -d"/" -f1 | rev | cut -d"." -f1)
# initialize creates the necessary directories and defines variables
/usr/bin/time -v bash "$basedir/scr/initialize.sh" >  "$basedir/scr/log/initialize.log" 2> "$basedir/scr/log/initialize.err"
# bwa indexing of the genome 
/usr/bin/time -v bash "$basedir/scr/indexing.sh" >  "$basedir/scr/log/indexing.log" 2> "$basedir/scr/log/indexing.err"

## 1) function to process individual sample.
process_yeast() {

source ./config

 local ind=$1
  
  if grep -qw "$ind" "$cpsdir/cps.txt"; then
    echo "Sample $ind already processed. skipping..."
    return 0
  fi
# bwa short read mapping
/usr/bin/time -v bash "$basedir/scr/mapping.sh" $ind >  "$basedir/scr/log/mapping.log" 2> "$basedir/scr/log/mapping.err"
# samtools parsing and genotyping 
/usr/bin/time -v bash "$basedir/scr/genotype.sh" $ind >  "$basedir/scr/log/genotype.log" 2> "$basedir/scr/log/genotype.err"
# handle filtered samples and positions 
/usr/bin/time -v bash "$basedir/scr/parsefilt.sh" $ind >  "$basedir/scr/log/parsefilt.log" 2> "$basedir/scr/log/parsefilt.err"
# get stats on gVCF for subsequnet assessment of the quality of the sample 
/usr/bin/time -v bash "$basedir/scr/getstats.sh" $ind >  "$basedir/scr/log/getstats.log" 2> "$basedir/scr/log/getstats.err"
# keep track of strains processed 
/usr/bin/time -v bash "$basedir/scr/cpcontroller.sh" $ind >  "$basedir/scr/log/cpcontroller.log" 2> "$basedir/scr/log/cpcontroller.err"

}
# export the function above 
export -f process_yeast
# run the function in parallel along the samples
find "$basedir/map" -name '*bam' | rev | cut -d"/" -f1 | rev | cut -d"." -f1 | sort -u | parallel -j "$nSamples" process_yeast {}
