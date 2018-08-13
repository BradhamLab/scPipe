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

### 1. Read Quality Control

`input`: raw reads (`.fastq.gz`)<br>
`output`: trimmed and filtered reads (`.fastq.gz`)

Perform quality control using the `fastp` tool [(link)](https://github.com/OpenGene/fastp) by trimming low quality regions and adapter sequences in reads, and filtering reads with too many ambiguous bases (Ns) or reads with low sequence complexity.

**bioRxiv Pre-Print**<br>
<sub>
Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu. fastp: an ultra-fast all-in-one FASTQ preprocessor. bioRxiv 274100; doi: https://doi.org/10.1101/274100
</sub>


### 2. Read Alignment

`input`: trimmed and filtered reads (`fastq.gz`)<br>
`output`: aligned reads (`.bam`, `.sam`)

Align filtered reads to the provided genome using `STAR`[(link)](https://github.com/alexdobin/STAR).

**Original Paper**<br>
<sub>
Dobin, A. Davis CA, Schlesinger F, Drenkow J. Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner.  Bioinformatics. 2013. 29. 1. pp 15-21.
</sub>


### 3. Expression Quantification via Read Counts

Retrieve fragment counts of paired end data using `featureCounts` [(link)](https://github.com/simon-anders/htseq). 

`input`: filtered alignments (`.bam`)<br>
`output`: raw read count matrix (`.csv`)

**Original Paper**<br>
<sub>
Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014. 30. 7. pp 923-930.
</sub>

### 4. Coverage Quality Control

Remove genes without any counts associated with it across cells. Remove cells with greater than 90% dropout. 

`input`: read counts (`.csv`)<br>
`ouput`: filtered counts (`.csv`)


### 5. Expression Matrix Normalization

`input`: raw read count matrix (`.csv`)<br>
`output`: within-sample normalized count matrix.

Normalize read counts using `SCnorm` [(link)](https://github.com/rhondabacher/SCnorm) if dropout is below 80%, otherwise use `scran` [(link)](http://bioconductor.org/packages/release/bioc/html/scran.html)

**Original Papers**<br>
<sub>
Bacher R, Chu LF, Leng N, Gasch AP, Thomson JA, Stewart RM, Newton M, Kendziorski C. SCnorm: robust normalization of single-cell RNA-seq data. Nature Methods. 2017 Jun 1;14(6):584-6.
</sub>

<sub>
Lun ATL, McCarthy DJ, Marioni JC. Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. Genome Biology. 2016 Feb 1;17(75).
</sub>

### 6. Batch Effect Removal

`input`: within-sample normalized count matrix.<br>
`output`: batched removed normalized count matrix.

Remove batch effects using mutual nearest neighbors [(link)](http://bioconductor.org/packages/release/bioc/html/scran.html).

**Orginal Paper**<br>
<sub>
Haghverdi L, Lun ATL, Mordan MD, Marioni JC. Batch effects in single-cell RNA-sequencing data are corrected by matching mutaul nearest neighbors. Nature Biotechnology. 2018. 26:421-427.
</sub>




<sup>*</sup> *This pipeline is currently being developed and does not exist in a complete/functional state.*
