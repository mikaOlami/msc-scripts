#!/bin/bash

# get chimeras links file, checks for valid chimeras
# link file is: Gene1 Start1 End1 Gene2 Start2 End2
#  USAGE: ~/scripts/ChimerasToLinks.sh [chimeras file]

chimeras_file=$1;

file=$(basename ${chimeras_file%.*})_links.out;

awk '{print $3 "\t" $4 "\t" $5 "\t" $10 "\t" $11 "\t" $12 "\t"}' $chimeras_file | awk '$2~/^[0-9]+/ && $3~/^[0-9]+/ $5~/^[0-9]+/ && $6~/^[0-9]+/ {print}' > $file;