#!/bin/bash

# 1. Enter the directory so relative paths work
cd /home/logan/software/v2p

# 2. Run with Singularity (-s)
#    -i: Relative path (file is in current dir)
#    -o: Relative path (output dir in current dir)
#    -a: Absolute path (required for the bind mount)
#    -s: TELLS THE SCRIPT TO USE SINGULARITY

bash get_predictions.sh \
  -i "pgpc4-clean.vcf" \
  -o "output2" \
  -c 8 \
  -a "$(pwd)/cadd_data/" \
  -s
