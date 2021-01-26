# AlternateReferenceMaker 
[![DOI](https://zenodo.org/badge/333140458.svg)](https://zenodo.org/badge/latestdoi/333140458)

is a command line based tool to integrate called variants into genomic sequences. \
To operate AlternateReferenceMaker needs the positional arguments: the FASTA file and the VCF, in this exact order. \
AlternateReferenceMaker can handle called variants of a multi-FASTA file.

## Requirements
[Python](https://www.python.org/downloads "Download Python") 2.7 or 3

## Usage
```bash
./AlternateReferenceMaker.py FASTA VCF 
```
### Positional Arguments:
```
  INPUT           name of input FASTA file containing FASTA sequences followed by name of VCF containing called variants
```
