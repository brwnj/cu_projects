#!/usr/bin/env python
# coding=utf-8
"""
"""

samples = ["MCF7-0", "MCF7-0_fraction", "MCF7-0_totalRNA", "MCF7-6",
            "MCF7-6_fraction", "MCF7-6_totalRNA", "PK12-0", "PK12-6"]
data_dir = "/vol1/home/brownj/projects/kabos_ribosome/data/common"
results_dir = "/vol1/home/brownj/projects/kabos_ribosome/results/common"

trim_adapter = "AGATCGGAAGAG"
trim_minlength = "25"

filter_idx = "/vol1/home/brownj/ref/hg19/rRNA"
filter_seedlen = 23

tophat_idx = "/vol1/home/brownj/ref/hg19/hg19"
tophat_segment_mismatches = 1
tophat_no_novel_juncs = True
tophat_transcriptome_index = "/vol1/home/brownj/ref/hg19/transcriptome/genes"
tophat_T = True

# rule all:
#     input: expand("/vol1/home/brownj/projects/kabos_ribosome/results/common/{sample}/Alignments/{sample}.bam".format(results=results_dir), sample=samples)

# rule filter:
#     input: rules.trim.output
#     output: "/vol1/home/brownj/projects/kabos_ribosome/data/common/$sample_trimmed.fastq.gz"
#
#     filt = ("{cat} {fastq} "
#             "| bowtie --un {result} {opts} - "
#             "> /dev/null "
#             "2> {stats}").format(cat="zcat" if sample.compressed else "cat",
#                                  fastq=sample.filename,
#                                  result=temp_file,
#                                  opts=opts,
#                                  stats=stats_file)
#     gzip = "gzip {fastq}".format(fastq=temp_file)
#     rename = "mv {temp}.gz {result}".format(temp=temp_file,
#                                             result=result_file)

rule trim:
    input: "/vol1/home/brownj/projects/kabos_ribosome/data/common/MCF7-0.fastq.gz"
    output: "/vol1/home/brownj/projects/kabos_ribosome/data/common/trimmed/{sample}.fastq.gz"
    params: TA="AGATCGGAAGAG", TML="25"
    shell: "zcat {input} | fastx_clipper -a {params.TA} -l {params.TML} -c -n -v -Q33 | fastx_trimmer -Q33 -f 2 | gzip -c > {output}"
