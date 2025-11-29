This repository contains instructions and code to run the Variant-to-Phenotype (V2P) variant effect predictor. For many users, accessing V2P 
via our website, v2p.ai, will be the best option. There, users may process up to 100,000 variants at a time. For advanced users or those who need 
higher throughput, offline installation may be a better option.

# Installing V2P

## Install git-lfs
https://git-lfs.com/

## Download the V2P repository

```
git clone https://github.com/davidfstein/v2p.git
```
Set the V2P working directory to the location that you cloned the repository.
```
export V2P_DIR=/path/to/dir
```
You can add this to your bash config if you would like it to persist.
```
echo "export V2P_DIR=/path/to/dir" >> ~/.bashrc
```

## Install dependencies

V2P's dependencies may be installed with conda. 
```
conda env create --file v2p.yaml
```
**Note: Some users report issues creating the conda environment. In some cases replacing the defaults channel with nodefault may fix this.

Activate the conda environment.
```conda activate v2p```
An additional python dependency must be installed after creating the conda environment.
This only needs to be done once. 
```
cd bin/scikit-multilearn
pip install .
```
V2P also requires an installation of docker or singularity to be available on the system. See here for installation instructions: https://docs.docker.com/engine/install/. If using docker, be sure to add your user to the docker group: https://docs.docker.com/engine/install/linux-postinstall/.

### Download V2P data

V2P relies on a large set of features to derive its predictions. These features are collected from Ensembl's VEP and other sources.
In total, downloading the feature data requires ~564GB of free disk space. 

Links to the V2P data will be emailed to you after registering on the v2p.ai website at https://www.v2p.ai/downloads

#### wget/curl download
Download the vep.tar.gz and hpo.db.gz files. 
```
cd v2p
# These should be placed in the V2P repository directory
tar xzf vep.tar.gz > ${V2P_DIR}/.vep
gunzip -c hpo.db.gz > ${V2P_DIR}/hpo.db
```

Download the cadddb.zip file.
```
# This can be placed anywhere on your system. You must provide the path to this data
# when running V2P
mkdir cadd_data
cd cadd_data
unzip cadddb.zip
mv cadddb/*db .
```

## (Optional) Download precomputed variants

Precomputed predictions for all posible single nucleotide variants and gnomAD indels are available for download. 
Installation will result in greatly increased speed for most variant sets. However, there is an additional 
disk space requirement.

#### wget/curl download
Download the snv_predictions.zip and indel_predictions.zip files.
```
# This can be placed anywhere on your system.
mkdir predictions
cd predictions 
unzip snv_predictions.zip
mv snv_predictions/*db .
unzip indel_predictions.zip
mv indel_predictions/*db .
```

# Testing installation
You can test your installation without downloading the full database files (the vep.tar.gz and hpo.db.gz must be downloaded).
```
cp ${V2P_DIR}/test/test.vcf ${V2P_DIR}/test.vcf
bash get_predictions.sh -i test.vcf -o results.csv -a ${V2P_DIR}/test/
python test/test.py results.csv
```
If the installation was successful, you will see "Generated predictions matched the expected predictions". 
(Runtime on CentOS 7.9.2009, 1 core: ~30s)
(Tested on CentOS 7.9.2009, Ubuntu 18.04, OSx)

# Running V2P

The input to V2P is a VCF file containing the variants you wish to score in hg38 coordinates. 
The VCF file must contain `#CHROM   POS ID  REF ALT QUAL    FILTER  INFO` in the header. 
Other header lines should not be included.

To run V2P provide the path to the VCF, a path where the output will be stored, and the path to the downloaded features.
You may optionally pass the path to the precomputed predictions, a number of CPUs to be used (default=1) and, 
if scoring variants from a single gene, you may provide a gene name to be prioritized by VEP. 
Run ```bash get_predictions.sh -h``` for a description of the available parameters.

Example
```
bash get_predictions.sh -i /path/to/input -o /path/to/output -a /absolute/path/to/annotations -p /path/to/precomputed_predictions -c 5
```
with singularity instead of docker
```
bash get_predictions.sh -i /path/to/input -o /path/to/output -a /absolute/path/to/annotations -p /path/to/precomputed_predictions -c 5 -s
```

# Alternate models

We also provide a version of V2P trained without other variant effect predictor outputs as features. 
You can run this version of the model using the -n flag. 

Example
```
bash get_predictions.sh -i /path/to/input -o /path/to/output -a /absolute/path/to/annotations -p /path/to/precomputed_predictions -c 5 -n
```
