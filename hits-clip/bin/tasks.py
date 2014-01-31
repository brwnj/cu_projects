#!/usr/bin/env python
# encoding: utf-8
"""
"""

from bsub import bsub
from invoke import task
from pathlib import Path
from string import Template
from yaml import dump, load

config_path = Path('/vol1/home/brownj/projects/hits-clip/bin/config.yaml')
assert config_path.exists(), "config file not found"
config = load(open(config_path.as_posix()))

fastqs = Path(config['fastqs'])
assert fastqs.is_dir()
results = Path(config['results'])
if not results.exists():
    results.mkdir(parents=True)

NOVOALIGN = Template(("novoalign -d $index -f $fastq -a -o SAM -r A 20 -e 100 -c $cpus -k 2> $summary "
                        "| samtools view -ShuF4 - "
                        "| samtools sort -o - $sample.temp -m 8G "
                        "> $bam"))

@task()
def trim():
    """trim reads by design using script that checks for sequence."""
    submit = bsub("trim", P=config['project_id'], verbose=True)
    for sample in config['trimmed_samples']:
        # name of the original file and the eventual result
        fastq = fastqs / "{sample}.fastq.gz".format(sample=sample)
        # rename
        untrimmed_fastq = "{base}.untrimmed.fastq.gz".format(base=fastq.absolute().as_posix().split(".fastq.gz")[0])
        fastq.rename(untrimmed_fastq)
        # trim
        python adapter_trim.py -a adapter -l 3 untrimmed_fastq | gzip -c > fastq
        python ../../bin/scripts/trim_adapter.py --reverse-complement -a CCGCTGGAAGTGACTGACAC PK12-24.untrimmed.fastq.gz | gzip -c > PK12-24.fastq.gz
        print cmd
        # submit(cmd)


@task()
def align():
    """align the reads using Novoalign."""
    cpus = str(config['align']['cpus'])
    novoidx = Path(config['align']['index'])
    assert novoidx.exists(), "index not found"

    submit = bsub("align", P=config['project_id'], n=cpus, R="select[mem>16] rusage[mem=16] span[hosts=1]", verbose=True)

    for sample in config['samples']:
        fastq = fastqs / '{sample}.fastq.gz'.format(sample=sample)

        # check results dir
        sample_results = results / sample
        if not sample_results.exists():
            sample_results.mkdir()

        bam = results / sample / '{sample}.bam'.format(sample=sample)
        if bam.exists(): continue

        alignment_summary = results / sample / '{sample}.alignment_stats.txt'.format(sample=sample)
        cmd = NOVOALIGN.substitute(index=novoidx.as_posix(),
                                    fastq=fastq.as_posix(),
                                    cpus=cpus,
                                    summary=alignment_summary,
                                    sample=sample,
                                    bam=bam)
        submit(cmd)


@task()
def removedups():
    # if [[ ! -f $rmdups_bam ]]; then
    #     samtools rmdup -s $bam $rmdups_bam
    # fi
    pass


@task()
def genomedata():
    """
    bams=$RESULTS/*/*bam
    for file in $bams; do
        if [[ ! -f $file.bai ]]; then
            samtools index $file
            bsub -J indexbam -o $file.out -e $file.err -P $PI -K "samtools index $file" &
        fi
    done
    wait

    if [[ ! -d $GENOMEDATA ]]; then
        # no need to wait for this to finish
        bam2gd.py $SIZES $FASTAS $bams -o $GENOMEDATA -p pillai_kabos_hitsclip &
    fi

    if the archive exists and there are tracks in there, need to add tracks
    one at a time!
    """
    pass
