#!/bin/bash

###########################################################
#			                 INITIALIZATION               #
#                                                         #
# This script initializes the environment for the sample  #
# processing pipeline. It loads configuration parameters, #
# sets up necessary directories, and checks dependencies. #
###########################################################

basedir=$(dirname "$(pwd)")

config_file="$basedir/scr/config"

set -a

source "$config_file"

set +a

# Create output directories
mkdir -p "$mapdir" "$tmpdir" "$logdir" "$cpsdir" "$statdir" "$gvcfdir"  "$bamstatdir" "$gvcfstatdir"

# Check dependencies
all_tools=("samtools" "bcftools" "bwa")

> "$logdir/tools"

> "$logdir/latest" 

echo "$(date)" > "$logdir/latest"

for tool in "${all_tools[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        echo "Error: $tool is not installed. Please install it before proceeding."
        exit 1
    else
    command -v "$tool" >> "$logdir/tools"
    fi
done

[ ! -f "$cpsdir/cps.txt" ] && touch "$cpsdir/cps.txt"

echo "initialization complete :)"
