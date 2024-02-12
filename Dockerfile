# Use the Miniforge3 image as a base, which comes with conda
FROM condaforge/miniforge3:latest

# Set the working directory in the container to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Create the conda environment from the YAML file
RUN mamba env create -f /app/envs/RNAseq.yaml

# Volumes for input data, configuration, and output results
VOLUME ["/data", "/results"]

# Create a reference directory
RUN mkdir -p /app/reference
RUN mkdir -p /app/annotation

# Download and decompress the reference genome and annotation
RUN wget -q -O /app/reference/homo_sapiens.GRCh38.dna.primary_assembly.fa.gz https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    && gunzip /app/reference/homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

RUN wget -q -O /app/annotation/Homo_sapiens.GRCh38.110.gtf.gz http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
    && gunzip /app/annotation/Homo_sapiens.GRCh38.110.gtf.gz

# Set an environment variable for the default number of cores
ENV SNAKEMAKE_CORES=1

# Ensure the container is executable
CMD ["/bin/bash"]
