# scPipe

Pipeline used for single-cell RNAseq read alignment in the Bradham Lab at Boston University.

The pipeline is implemented using SnakeMake.

## Conda Environments

Two conda environments are required to run the pipeline:

1. Alignment

Create the alignment environment using the `requirements.txt`. In a terminal, with an accessible `conda` installation, issue the following command:

```{bash}
conda create --name alignment --file alignment_spec.txt 
```

2. MultiQC

The MultiQC environment should be built using the following instructions:

```{bash}
conda create -n multiqc --no-default-packages
conda install pip
pip install --upgrade --force-reinstall git+https://github.com/ewels/MultiQC.git --ignore-installed certifi
```
