# scPipe

Pipeline used for single-cell RNAseq read alignment in the Bradham Lab.

The pipeline is implemented using SnakeMake.

## Conda Environments

Two conda environments are required to run the pipeline

1. Alignment

2. MultiQC

The MultiQC environment should be built using the following instructions

```{bash}
conda create -n multiqc --no-default-packages
conda install pip
pip install --upgrade --force-reinstall git+https://github.com/ewels/MultiQC.git --ignore-installed certifi
```
