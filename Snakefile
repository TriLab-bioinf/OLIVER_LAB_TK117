# vim: set ft=python:

# RNAseq workflow v1.0
# Hernan Lorenzi
# hernan.lorenzi@nih.gov
# Workflow requires to have a STAR index already available within the data/00ref directory
# Also, it is necessary to configure the config.yml file accordingly to include all metadatata required.

import os
import glob

configfile: "config/config.yaml"
samples = config["samples"].keys()
genome = config["reference"]["genome_file"]
annotation = config["reference"]["ensembl_gtf"]

# Build directory structure
paths = ["data/00ref", "data/00reads", "data/00adapters",
        "results/01trim", "results/02abundant","results/03map_reads",
        "results/03map_reads2","results/04dedup", "results/05bigwig",
        "results/06fastqc_raw",
        "results/06fastqc_trim","results/07multiqc"
        ]
for path in paths:
    os.makedirs(path, exist_ok = True)

# Functions
def get_fq1(wildcards):
            return [ my_files for my_files in glob.glob(f"data/00reads/{wildcards.sample}.R1.fastq.gz")]

def get_fq2(wildcards):
            b = [ my_files for my_files in glob.glob(f"data/00reads/{wildcards.sample}.R1.fastq.gz")]
            c = list(map(lambda x: str.replace(x, ".R1", ".R2"), b))
            return c

# Set what rules to run locally
localrules: all #,
            #build_abundant_db

rule all:
    # IMPORTANT: output file fo all rule has to match the name specified in the output file
    # and not include suffixes that the command use might add to it.
    input:  #o10 = expand("results/01trim/{s}.{e}.fastq.gz", e=["1P","2P","1U","2U"], s=samples)
            #o0 = "results/05counts/read_counts.summary",
            o1 = expand("results/06fastqc_raw/{s}.R1_fastqc.html", s=samples),
            o2 = expand("results/06fastqc_raw/{s}.R2_fastqc.html", s=samples),
            o3 = expand("results/06fastqc_trim/{s}.1P_fastqc.html", s=samples),
            o4 = expand("results/06fastqc_trim/{s}.2P_fastqc.html", s=samples),
            o7 = "results/07multiqc/multiqc_done.flag",
            o8 = expand("results/03map_reads2/{s}.sam", s=samples),
            #o8 = "results/05correlation/multiBamSummary.results.npz",
            o9 = expand("results/05bigwig/{s}.bw", s=samples)
            #o11 = "data/00ref/SA"
            #o12 = expand("results/01trim/{s}.U.fastq.gz", s=samples),
            #o13 = expand("results/00merged_reads/{s}.R1.fastq.gz", s=samples),
            #o14 = expand("results/00merged_reads/{s}.R2.fastq.gz", s=samples),
            #expand("00map_reads/{s}.", s=samples) #,
            #expand("00abundant/{s}.fastq.1.gz", s=samples),
            #expand("00abundant/{s}.fastq.2.gz", s=samples)

rule merge_rep:
    input: fq1 = get_fq1,
           fq2 = get_fq2
    output: merge1 = temp("results/00merged_reads/{sample}.R1.fastq.gz"),
            merge2 = temp("results/00merged_reads/{sample}.R2.fastq.gz")
    resources: 
        cpu_per_task = 1,
        partition = "quick",
        time = "3:00:00"
    threads: 1
    shell:
        """
        cat {input.fq1}  > {output.merge1}
        cat {input.fq2}  > {output.merge2}
        """

rule trimming:
    input:  fq1 = "results/00merged_reads/{sample}.R1.fastq.gz",
            fq2 = "results/00merged_reads/{sample}.R2.fastq.gz", 

    output:
            fq1P = "results/01trim/{sample}.1P.fastq.gz",
            fq2P = "results/01trim/{sample}.2P.fastq.gz",
            fqU = "results/01trim/{sample}.U.fastq.gz"
    params: 
            "ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20 overwrite=t"
    resources:
        cpus_per_task = 16,
        partition = "quick",
        time = "4:00:00"
    threads: 16
    log:    log1 = "results/01trim/{sample}.log",
            log2 = "results/01trim/{sample}.stats.log"
    benchmark:
            "benchmarks/trim/{sample}.tsv"
    shell:
        """
        bbduk.sh -Xmx1g threads={threads} \
            in1={input.fq1} in2={input.fq2} \
            out1={output.fq1P} out2={output.fq2P} outs={output.fqU} \
            ref=data/00adapters/truseq.fa.gz \
            {params} stats={log.log2} 2> {log.log1}  
        """


       
rule map_reads:
    input: fq1 = "results/01trim/{sample}.1P.fastq.gz", #"results/02abundant/{sample}.fastq.1.gz",
           fq2 = "results/01trim/{sample}.2P.fastq.gz" #"results/02abundant/{sample}.fastq.2.gz",
    output: "results/03map_reads/{sample}.bam"
    threads: 16
    resources:
        cpus_per_task = 16,
        partition = "norm",
        time = "14:00:00",
        mem_mb = 32000,
        gres = "lscratch:20"
    benchmark:
        "benchmarks/map_reads/{sample}.tsv"
    params: genome = "data/00ref/Dmel.BDGP6.46.dna_rm.toplevel"
    shell:
        """
        bowtie2 -p 16 \
        -x {params.genome} \
        -1 {input.fq1} \
        -2 {input.fq2} \
        | samtools view -@8 -Sb -|samtools sort -@8 - -o {output}
        """


rule map_reads2:
    input: fq1 = "results/01trim/{sample}.1P.fastq.gz", #"results/02abundant/{sample}.fastq.1.gz",
           fq2 = "results/01trim/{sample}.2P.fastq.gz" #"results/02abundant/{sample}.fastq.2.gz",
    output: "results/03map_reads2/{sample}.sam"
    threads: 16
    resources:
        cpus_per_task = 16,
        partition = "norm",
        time = "14:00:00",
        mem_mb = 32000,
        gres = "lscratch:20"
    benchmark:
        "benchmarks/map_reads2/{sample}.tsv"
    params: genome = "data/00ref/Dmel.BDGP6.46.dna_rm.toplevel"
    shell:
        """
        bwa-mem2 mem -t 16 \
        {params.genome} \
        {input.fq1} \
        {input.fq2} \
        -o {output}
        """

rule remove_duplicates:
    input: "results/03map_reads/{sample}.bam"
    output: "results/04dedup/{sample}.sorted.dedup.bam"
    params: "READ_NAME_REGEX=null REMOVE_DUPLICATES=false"
    log: "results/04dedup/{sample}.sorted.dedup.metrics.txt"
    threads: 24
    benchmark:
        "benchmarks/remove_duplicates/{sample}.tsv"
    resources:
        cpus_per_task = 4,
        mem_mb = 96000,
        partition = "quick",
        time = "4:00:00",
        gres = "lscratch:20"
    shell:
        """
        picard -Xmx16g -Xms16g -XX:ParallelGCThreads=5 MarkDuplicates \
         I={input} \
         O={output} \
         M={log} \
         {params}
        samtools index {output}
        """

rule make_bigwig:
    input: "results/04dedup/{sample}.sorted.dedup.bam"
    output: "results/05bigwig/{sample}.bw"
    params: "--binSize 10 --normalizeUsing BPM" #  + "--filterRNAstrand [forward/reverse]" to plot strand-specific data
    threads: 8
    resources:
        cpus_per_task = 8,
        partition = "norm",
        time = "14:00:00"
    shell:
        """
        bamCoverage -p {threads} -b {input} -o {output} {params}
        """

rule fastqc:
    input: raw1 = "results/00merged_reads/{sample}.R1.fastq.gz",
           raw2 = "results/00merged_reads/{sample}.R2.fastq.gz",
           trim1p = "results/01trim/{sample}.1P.fastq.gz",
           trim_u = "results/01trim/{sample}.U.fastq.gz",
           trim2p = "results/01trim/{sample}.2P.fastq.gz"
    output: 
            o1 = "results/06fastqc_raw/{sample}.R1_fastqc.html",
            o2 = "results/06fastqc_raw/{sample}.R2_fastqc.html",
            o3 = "results/06fastqc_trim/{sample}.1P_fastqc.html",
            o5 = "results/06fastqc_trim/{sample}.2P_fastqc.html"
    threads: 8
    resources:
        cpus_per_task = 8,
        partition = "quick",
        time = "4:00:00",
        mem_mb = 4000
    params:
        "--quiet"
    shell:
        """
        fastqc {params} -t {threads} -o results/06fastqc_raw {input.raw1} {input.raw2}
        fastqc {params} -t {threads} -o results/06fastqc_trim {input.trim1p} {input.trim2p}
        """
absolute_path = "/gpfs/gsfs12/users/wangy80/TK117/workflow/"

rule multiqc:
    input: 
           i1 = "results/01trim",
           i2 = "results/03map_reads", 
           i3 = "results/03map_reads2",
           i4 = "results/04dedup",
           i7 = "results/06fastqc_raw",
           i8 = "results/06fastqc_trim"
    output: "results/07multiqc/multiqc_done.flag"
    resources:
        partition = "quick",
        time = "4:00:00",
        mem_mb = 4000                    
    shell:
        """
        multiqc -f -d -o results/07multiqc {input.i1} \
                {input.i7} \
                {input.i8} \
                {input.i3} \
                {input.i4} 
        touch {output}
        """