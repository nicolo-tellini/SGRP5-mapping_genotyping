#!/bin/bash
###########################################################
#                        MAPPING                          #
#                                                         #
# This script performs the mapping of the short reads and #
# the postprocessing of the BAM file.                     #
###########################################################

 source ./config

 ind=$1

 forward_read=$(find "$seqdir" -name "${ind}_1."*"fastq.gz" -print -quit)

 reverse_read=$(find "$seqdir" -name "${ind}_2"*"fastq.gz" -print -quit)
 
 sam_file="$mapdir/${ind}.sam"
 
 genome=$(find "$(realpath $basedir/rep)" -name "*.fa")
 
 ref=$(echo $genome | rev | cut -d"/" -f1 | rev | cut -d"." -f1)
 
 bam_file="$mapdir/${ind}.${ref}.bam"

 # short-read alignment
 if [ -n "$reverse_read" ]; then
	 
	 bwa mem -t "$nThreads" $genome "$forward_read" "$reverse_read" -o "$sam_file" 2> /dev/null
 else
	 bwa mem -t "$nThreads" $genome "$forward_read" -o "$sam_file" 2> /dev/null
 fi

 
 ## SAM processing
 samtools fixmate -@$nThreads -O bam,level=1 -m "$sam_file" "$mapdir/${ind}.fix.bam"
 rm -f "$sam_file"
  
 samtools sort -l 1 -@$nThreads "$mapdir/${ind}.fix.bam" -T "$tmpdir/${ind}.tmp.bam" -o "$mapdir/${ind}.fix.srt.bam"
 rm -f "$mapdir/${ind}.fix.bam"
  
 samtools markdup -@$nThreads -O bam,level=1 "$mapdir/${ind}.fix.srt.bam" "$mapdir/${ind}.fix.srt.mrk.bam"
 rm -f "$mapdir/${ind}.fix.srt.bam"
  
 samtools view -@$nThreads "$mapdir/${ind}.fix.srt.mrk.bam" -o "$bam_file"
 rm -f "$mapdir/${ind}.fix.srt.mrk.bam"
  
 samtools index -@2 "$bam_file"
