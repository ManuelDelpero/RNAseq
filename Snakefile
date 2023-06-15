configfile: "config.yaml"

SAMPLES = config['samples']
output_dir = config['output_dir']
raw_dir = config['rawdir']
alignment_dir = output_dir + '/alignment'
log_dir = output_dir + '/logs'
benchmark_dir = log_dir + '/benchmarks'
alignment_dir = output_dir + '/alignment'
qc_dir = output_dir + '/qc'
stats_dir = output_dir + '/stats'
gene_annotation_file = config['gtf_annotation_file']
ref_fasta = config['ref']
ref_dir = output_dir + '/reference'
count_dir = output_dir + '/counts'
num_threads = config['computing_threads']

results_counts = expand(count_dir + '/{sample}_count.txt', sample = SAMPLES)
results_fastqc = expand(qc_dir + '/fastqc/{sample}_fastqc.html', sample = SAMPLES)

rule all:
    input:
        results_counts + 
        results_fastqc

rule fastp:
    """
    Trim reads and perform QC using fastp with raw fastq
    """
    input:
        R1 = raw_dir + '/{sample}_R1.fastq.gz',
        R2 = raw_dir + '/{sample}_R2.fastq.gz'
    output:
        R1 = qc_dir + '/fastp/{sample}_R1_fastp.fastq.gz',
        R2 = qc_dir + '/fastp/{sample}_R2_fastp.fastq.gz',
        stats = qc_dir + '/fastp/{sample}_fastp.html'
    log:
        log_dir + '/{sample}.fastp.log'
    threads:
	    num_threads
    benchmark:
        benchmark_dir + '/fastp/{sample}.tsv'
    shell:
        'fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} -h {output.stats} --thread {threads} &> {log}'


rule hisat2_index:
    """
    index for Hisat2
    """
    input:
        config['ref']
    output:
        ref_dir
    log:
        log_dir + '/Hisat2_index.log'
    threads:
	    num_threads
    benchmark:
        benchmark_dir + '/Hisat2/Hisat2_index.tsv'
    shell:
        'hisat2-build {input} {output} -p {threads} &> {log}'


rule hisat2:
    """
    mapping and alignment using Hisat2
    """
    input:
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
        ref = ref_dir
    output:
        alignment_dir + '/hisat2/{sample}.sam'
    log:
        log_dir + '/{sample}.Hisat2.log'
    threads:
	    num_threads
    benchmark:
        benchmark_dir + '/Hisat2/{sample}.tsv'
    shell:
        """
        hisat2 -x {input.ref} -1 {input.R1} -2 {input.R2} -S {output} -p {threads} &> {log}
        """
		
rule sort:
    """
    Conversion and position sorting
    """
    input:
        rules.hisat2.output
    output:
        alignment_dir + '/hisat2/{sample}_sorted.bam'
    log:
        log_dir + '/{sample}.samtools_sort.log'
    threads:
	    num_threads
    benchmark:
        benchmark_dir + '/samtools_sort/{sample}.tsv'
    shell:
        """
        samtools sort {input} -o {output} -@ {threads} -m 16G &> {log}
        """

rule fastqc:
    """
    Run QC for each fastq
    """
    input:
        rules.sort.output
    output:
        qc_dir + '/fastqc/{sample}_fastqc.html'
    log:
        log_dir + '/{sample}_fastqc.log'
    threads:
	    num_threads
    benchmark:
        benchmark_dir + '/fastqc/{sample}.tsv'
    shell:
        'fastqc -o $(dirname {output}) {input} &> {log}'
		
rule Count:
    """
    Calculate raw reads count with HTseq
    """
    input:
        bam = rules.sort.output,
        annotation = gene_annotation_file
    output:
        count_dir + '/{sample}_count.txt'
    benchmark:
        benchmark_dir + '/HTseq/{sample}.tsv'
    shell:
        """
        htseq-count -f bam -r pos -t exon -i gene_id -s no {input.bam} {input.annotation} &> {output}
        """