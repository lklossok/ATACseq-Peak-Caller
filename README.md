# ATACseq-Peak-Caller
ATAC-seq processing pipeline implemented using Snakemake. The workflow downloads paired-end FASTQ files from the European Nucleotide Archive, performs read quality control and trimming, aligns reads to a user-specified reference genome, removes PCR duplicates, and generates genome-wide chromatin accessibility peaks.

# Pipeline steps:

download_fastq - Downloads paired-end FASTQ files for each SRR accession directly from ENA using wget, storing them under data/. File pairing is inferred from the URL list.

quality_control - Performs quality control trimming using fastp. Produces cleaned R1/R2 FASTQ files along with HTML and JSON QC reports.

bwa_mem - Aligns trimmed reads to the reference genome using BWA-MEM, pipes the SAM output into samtools sort, and generates a sorted, indexed BAM file for each sample.

deduplicate - Removes PCR duplicates using a multi-step samtools workflow. Produces a cleaned, duplicate-free BAM file.

index_dedup - Creates a BAM index (.bai) for each deduplicated BAM using samtools index. 

peak_calling - Calls ATAC-seq peaks with MACS2 using ATAC-specific parameters and outputs narrowPeak files for each accession.

# Inputs:

paired-end ATACseq FASTQ files from the European Nucleotide Archive 
hg38 bwa inex

# Outputs:

Processed ATACseq BAM files
Genome-wide narrowPeak files
