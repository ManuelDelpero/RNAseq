# Scalable RNAseq pipeline in Snakemake

This Snakemake pipeline processes and analyzes RNA sequencing data using various tools such as `fastp`, `hisat2`, `samtools`, `fastqc`, and `htseq-count`. The pipeline is configured using the `config.yaml` file, which contains the necessary parameters for running the pipeline. If more than 3 samples per group are provided the pipeline will also automatically run differential gene expression analysis using DESeq2.

## Prerequisites

Make sure you have Mamba installed on your system. If you don't have Mamba installed, you can follow the instructions at [mamba.io](https://mamba.readthedocs.io/en/latest/installation.html) to install Mamba.

## Pipeline Setup

1. Clone or download the pipeline code from the repository.
2. Ensure that the `config.yaml` file is present in the main directory of the repository, along with the Snakefile.

## Creating the Conda Environment

1. Create a new environment named "RNAseq" and install the dependencies from the `environment.yaml` file by running the following command:
   ```bash
   mamba env create -f envs/RNAseq.yaml
   
This will create a new environment called "RNAseq" and install all the required dependencies.

## Activating the Conda Environment and run the pipeline

   ```bash
   conda activate RNAseq
   snakemake --cores <num_cores>


## Output Directory Structure

The pipeline generates the following output directory structure:

output/
├── alignment/
│   └── hisat2/
│       ├── simulated1.bam
│       ├── simulated1.bam.bai
│       ├── simulated1.sam
│       ├── simulated2.bam
│       ├── simulated2.bam.bai
│       └── simulated2.sam
├── counts/
│   ├── simulated1_count.txt
│   └── simulated2_count.txt
├── qc/
│   ├── fastp/
│   │   ├── simulated1_R1_fastp.fastq.gz
│   │   ├── simulated1_R2_fastp.fastq.gz
│   │   ├── simulated2_R1_fastp.fastq.gz
│   │   └── simulated2_R2_fastp.fastq.gz
│   ├── fastqc/
│   │   ├── simulated1_fastqc.html
│   │   └── simulated2_fastqc.html
│   └── qualimap/
│       ├── simulated1_qualimap.html
│       └── simulated2_qualimap.html
├── logs/
│   ├── simulated1.fastp.log
│   ├── simulated1.Hisat2.log
│   ├── simulated1.samtools_sort.log
│   ├── simulated1_fastqc.log
│   ├── simulated2.fastp.log
│   ├── simulated2.Hisat2.log
│   ├── simulated2.samtools_sort.log
│   └── simulated2_fastqc.log
├── benchmarks/
│   ├── fastp/
│   │   ├── simulated1.tsv
│   │   └── simulated2.tsv
│   ├── Hisat2/
│   │   ├── Hisat2_index.tsv
│   │   ├── simulated1.tsv
│   │   └── simulated2.tsv
│   └── samtools_sort/
├── dex/
│   ├── DE_analysis.csv
│   ├── Normalized_counts.csv
│   ├── Heatmap.pdf
│   └── Volcano_Plot.pdf
└── multiqc/
    └── multiqc_report.html
