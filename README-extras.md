### Notes

This fork contains a few useful extras and was tested on a 24 CPU i9 + Blackwell RTX6000 / Ubuntu 24.04 system. I was unable to get the docker version to work. However, the singularity version worked fine.

**deepvariant/vcf_parser.py** -- converts either a gzipped or regular .vcf file into a new vcf file that is compatible with v2p by removing any vcf entries from unsupported chromosomes and random contigs that may interfere with processing, stripping unused fields, remove low quality + RefCall entries, and converting the CHROM entry to just a number or [1-22,X,Y]. This script was used to clean up a .vcf file from a genome that I processed with the DeepVariant module within the _nf-core/sarek_ suite. 

**runme_singularity.sh**  -- a shell script that starts the analysis

**runme_singularity_predictions.sh** -- a shell script that uses the precomputed databases

**scripts/getPrecomputed.py** -- modified the original script to provide some useful debugging information during processing


