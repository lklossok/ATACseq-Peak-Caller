import os
import re

configfile: "config/config.yaml"

# Read URLs from file
with open(config["ENA_dir"]) as f:
    URLS = [line.strip() for line in f if line.strip()]

# Extract pairs using regex SRRxxxx_{1,2}.fastq.gz
pairs = {}

for url in URLS:
    fname = os.path.basename(url)

    m = re.match(r"(SRR\d+)_(1|2)\.fastq\.gz", fname)
    if not m:
        raise ValueError(f"Unexpected filename format: {fname}")

    accession, mate = m.groups()
    pairs.setdefault(accession, {})[f"r{mate}"] = {
        "url": url,
        "file": fname
    }

ACCESSION = sorted(pairs.keys())

INDEX = config["index_dir"]
BWA_INDEX = f"{INDEX}/genome.fa"


rule all:
    input:
        expand("data/narrowPeak/{accession}_peaks.narrowPeak", accession=ACCESSION)


rule download_fastq:
    output:
        temp("data/{accession}_R1.fastq.gz"),
        temp("data/{accession}_R2.fastq.gz")
    run:
        urls = pairs[wildcards.accession]
        shell(f"wget -O {output[0]} {urls['r1']['url']}")
        shell(f"wget -O {output[1]} {urls['r2']['url']}")

rule quality_control: 
    input:
        r1 = "data/{accession}_R1.fastq.gz",
        r2 = "data/{accession}_R2.fastq.gz"
    output:
        r1 = temp("data/{accession}_R1.trimmed.fastq.gz"),
        r2 = temp("data/{accession}_R2.trimmed.fastq.gz"),
        html = "data/fastp/{accession}.fastp.html",
        json = "data/fastp/{accession}.fastp.json"
    threads: 8
    conda: "envs/ATACseq.yaml"
    shell:
        """
        mkdir -p data/fastp

        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            -h {output.html} \
            -j {output.json} \
            --thread {threads} \
            --qualified_quality_phred 30
        """

rule bwa_mem:
    input:
        r1= "data/{accession}_R1.trimmed.fastq.gz",
        r2= "data/{accession}_R2.trimmed.fastq.gz",
        idx= expand(f"{BWA_INDEX}.{{ext}}", ext=["amb","ann","bwt","pac","sa"])
    output:
        bam= temp("data/{accession}.sorted.bam")
    threads: 8
    conda: "envs/ATACseq.yaml"
    shell:
        """
        bwa mem -t {threads} {BWA_INDEX} {input.r1} {input.r2} |
        samtools sort -o {output.bam}
        samtools index {output.bam}
        """

rule deduplicate:
    input:
        bam = "data/{accession}.sorted.bam"
    output:
        dedup = "data/bam/{accession}.dedup.bam"
    conda: "envs/ATACseq.yaml"
    threads: 8
    shell:
        """
        mkdir -p data/bam

        samtools sort -n -@ {threads} -o {wildcards.accession}.namesort.bam {input.bam}
        samtools fixmate -m -@ {threads} {wildcards.accession}.namesort.bam {wildcards.accession}.fixmate.bam
        samtools sort -@ {threads} -o {wildcards.accession}.positionsort.bam {wildcards.accession}.fixmate.bam
        samtools markdup -@ {threads} -r {wildcards.accession}.positionsort.bam {output.dedup}
        rm {wildcards.accession}.namesort.bam {wildcards.accession}.fixmate.bam {wildcards.accession}.positionsort.bam
        """

rule index_dedup:
    input:
        "data/bam/{accession}.dedup.bam"
    output:
        "data/bam/{accession}.dedup.bam.bai"
    conda: "envs/ATACseq.yaml"
    shell:
        """
        samtools index {input}
        """

rule peak_calling:
    input:
        bam = "data/bam/{accession}.dedup.bam"
    output:
        narrowpeak = "data/narrowPeak/{accession}_peaks.narrowPeak"
    conda: "envs/ATACseq.yaml"
    shell:
        """
        mkdir -p data/narrowPeak

        macs2 callpeak -t {input.bam} \
            -f BAM \
            -n {wildcards.accession} \
            --outdir data/narrowPeak/ \
            -g hs \
            --nomodel \
            --shift -100 \
            --extsize 200 \
            -q 0.01
        """
