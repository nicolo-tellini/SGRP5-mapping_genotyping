#!/bin/bash

###########################################################
#                        INDEXING                         #
#                                                         #
# This script performs indexing  required for the         #
# mapping. It indexes the reference genome                #
###########################################################

source ./config

cd "$basedir/rep"

bwa index *.fa
