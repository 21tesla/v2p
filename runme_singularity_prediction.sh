#!/bin/bash

# 1. Enter the directory so relative paths work and we find the scripts
cd /home/logan/software/v2p

# 2. Run with Precomputed variants (-p)
#    -p: Must be a path valid on YOUR machine (Host). 
#        $(pwd)/predictions/ expands to /home/logan/software/v2p/predictions/

bash get_predictions.sh \
  -i "sample_full.vcf" \
  -o "output1" \
  -c 16 \
  -a "$(pwd)/cadd_data/" \
  -p "$(pwd)/predictions/" \
  -s
