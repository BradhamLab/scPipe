# scPipe

Pipeline used for single-cell RNAseq read alignment in the Bradham Lab at Boston University.

The pipeline is implemented using SnakeMake<sup>*</sup>.

## Conda Environments

Two conda environments are required to run the pipeline:

1. Alignment

Create the alignment environment using the `requirements.txt`. In a terminal, with an accessible `conda` installation, issue the following command:

```{bash}
conda create -n alignment --file alignment_spec.txt 
```

2. MultiQC

The MultiQC environment should be built using the following instructions<sup>*</sup>:

```{bash}
conda create -n multiqc pip --no-default-packages
source activate multiqc
pip install --upgrade --force-reinstall git+https://github.com/ewels/MultiQC.git --ignore-installed certifi
```

You may also need to install `Cython` for some package use in the `multiqc` environment.

This can be installed using the `conda install cython` command.

<sup>*</sup>Note, on some systems the Python 3 version of `MultiQC` fails due to the `click` library failing to deal with strings properly. If this is the case, specify `python=2.7` upon environment creation.


## Pipeline

This pipeline performs the necessary operations to take single-cell RNAseq data from paired-end raw reads to a normalized expression matrix. This transformation is done using the following tools/steps.

#### 1. Read Quality Control

input: raw reads (`.fastq`/`.fastq.gz`)<br>
output: trimmed reads (`.fastq`/`.fastq.gz`)

Raw reads are checked for high quality

#### 2. Read Alignment


#### 3. Alignment Quality Control


#### 4. Expression Quantification with Mapped Reads


#### 5. Expression Matrix Normalization

<sup>*</sup> *This pipeline is currently being developed and does not exist in a complete/functional state.*
