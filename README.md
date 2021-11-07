# USAGE

## Requirements

- conda 4.10.3

## Install environment

conda env create -f requirements.yml

conda activate CS5824

## Update environment

conda env export > requirements.yml

# Source
Full viral genomes have been download in Nov 7, 2021 from https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/

- viral.1.1.genomic.fna.gz
- viral.2.1.genomic.fna.gz
- viral.3.1.genomic.fna.gz
- viral.4.1.genomic.fna.gz

All is concatenated into a single ***full_viral.fna*** file with ***annotation.txt*** file containing the accession number and corresonding description.