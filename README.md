# Scalable RNAseq pipeline in Snakemake

This Snakemake pipeline processes and analyzes RNA sequencing data using various tools such as `fastp`, `hisat2`, `samtools`, `fastqc`, and `htseq-count`. The pipeline is configured using the `config.yaml` file, which contains the necessary parameters for running the pipeline.

## Prerequisites

Make sure you have Conda installed on your system. If you don't have Conda installed, you can follow the instructions at [conda.io](https://conda.io/projects/conda/en/latest/user-guide/install/) to install Conda.

## Pipeline Setup

1. Clone or download the pipeline code from the repository.
2. Ensure that the `config.yaml` file is present in the main directory of the repository, along with the Snakefile.

## Creating the Conda Environment

1. Install Mamba, a faster alternative to Conda, by running the following command:
   ```bash
   conda install mamba -n base -c conda-forge
2. Create a new Conda environment named "RNAseq" and install the dependencies from the `environment.yaml` file by running the following command:
   ```bash
   mamba env create -f envs/RNAseq.yaml
   
This will create a new environment called "RNAseq" and install all the required dependencies.

## Activating the Conda Environment

Before running the pipeline, you need to activate the "RNAseq" environment. You can do this by running the following command:

   ```bash
   conda activate RNAseq

## Running the Pipeline

To run the pipeline, execute the following command:
    ```bash
    snakemake --cores <num_cores>
