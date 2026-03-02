#!/bin/bash

set -x 

source ./config

genome=$(find "$(realpath $basedir/rep)" -name "*.fa")

ref=$(echo $genome | rev | cut -d"/" -f1 | rev | cut -d"." -f1)

/usr/bin/time -v bash "$basedir/scr/initialize.sh" >  "$basedir/scr/log/initialize.log" 2> "$basedir/scr/log/initialize.err"

/usr/bin/time -v bash "$basedir/scr/indexing.sh" >  "$basedir/scr/log/indexing.log" 2> "$basedir/scr/log/indexing.err"

## 1) function to process individual sample.
process_yeast() {

source ./config

 local ind=$1
  
  if grep -qw "$ind" "$cpsdir/cps.txt"; then
    echo "Sample $ind already processed. skipping..."
    return 0
  fi

/usr/bin/time -v bash "$basedir/scr/mapping.sh" $ind >  "$basedir/scr/log/mapping.log" 2> "$basedir/scr/log/mapping.err"

/usr/bin/time -v bash "$basedir/scr/genotype.sh" $ind >  "$basedir/scr/log/genotype.log" 2> "$basedir/scr/log/genotype.err"

/usr/bin/time -v bash "$basedir/scr/parsefilt.sh" $ind >  "$basedir/scr/log/parsefilt.log" 2> "$basedir/scr/log/parsefilt.err"

/usr/bin/time -v bash "$basedir/scr/getstats.sh" $ind >  "$basedir/scr/log/getstats.log" 2> "$basedir/scr/log/getstats.err"

/usr/bin/time -v bash "$basedir/scr/cpcontroller.sh" $ind >  "$basedir/scr/log/cpcontroller.log" 2> "$basedir/scr/log/cpcontroller.err"

}

export -f process_yeast

find "$basedir/map" -name '*bam' | rev | cut -d"/" -f1 | rev | cut -d"." -f1 | sort -u | parallel -j "$nSamples" process_yeast {}
