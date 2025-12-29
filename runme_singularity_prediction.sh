#!/bin/bash

# 1. Enter the directory so relative paths work and we find the scripts
cd /home/logan/software/v2p

# 2. Run with Precomputed variants (-p)
#    -p: Must be a path valid on YOUR machine (Host). 
#        $(pwd)/predictions/ expanded to /home/logan/software/v2p/predictions/ in my case

bash get_predictions.sh \
  -i "pgpc4-clean.vcf" \
  -o "pgpc4-output.parquet" \
  -c 16 \
  -a "$(pwd)/cadd_data/" \
  -p "$(pwd)/predictions/" \
  -s
