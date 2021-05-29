# scSimClassify

This repository contains a Python implementation of scSimClassify.

Four major modules are included:

- k-mer matrix generation module.

- CKG generation module.

- Cell matrix generation module.

- Classification module.

## Prerequisite

python 3.6

python packages: simhash, keras, sklearn, pandas, numpy, Bio


## Major module function

### k-mer matrix generation module
Both BAM and FASTA file can be the inputs. Sparse k-mer matrix is build by several steps, such as indexing k-mer, convert k-mer sequence to k-mer id, and generate k-mer matrix from jellyfish files.

### CKG generation module
By selecting informative k-mers and simhashing k-mer vectores into fingerprints, a CKG look-up table is generated.

### Cell matrix generation module
Based on CKG look-up table, the features of cells is transformed to the expression of CKGs from abundance of k-mers.

### Classification module
Four classification methods are provided in this module.


## Dataset description
Cell labels and examples of gene expression features in the experiments are provided.

## Usage instruction
Please modify the path to use the modules.








