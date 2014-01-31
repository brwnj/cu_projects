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

TRIM = Template(("python {script} --reverse-complement -a {adapter} "
                    "$untrimmed | gzip -c > $trimmed").format(\
                        script=config['trim_script'],
                        adapter=config['adapter_sequence']))

NOVOALIGN = Template(("novoalign -d $index -f $fastq -a -o SAM -r A 20 -e 100 -c $cpus -k 2> $summary "
                        "| samtools view -ShuF4 - "
                        "| samtools sort -o - $sample.temp -m 8G "
                        "> $bam"))

@task()
def trim_sequences():
    """
    trim reads by design using script that checks for sequence.

    writes new fastq will same name as previous
    renames untrimmed fastq $sample.untrimmed.fastq.gz
    """
    submit = bsub("trim", P=config['project_id'], verbose=True)
    for sample in config['trimmed_samples']:
        # name of the original file and the eventual result
        fastq = fastqs / "{sample}.fastq.gz".format(sample=sample)
        # renamed fastq
        untrimmed_fastq = Path("{base}.untrimmed.fastq.gz".format(base=fastq.as_posix().split(".fastq.gz")[0]))
        # check if complete
        if untrimmed_fastq.exists() and fastq.exists(): continue

        fastq.rename(untrimmed_fastq.as_posix())
        # trim
        cmd = TRIM.substitute(untrimmed=untrimmed_fastq.as_posix(), trimmed=fastq.as_posix())
        submit(cmd)


@task()
def align():
    """
    align the reads using Novoalign.

    writes $sample.bam to $results/$sample/$sample.bam
    also writes an alignment summary in the same dir as
    $sample.alignment_stats.txt
    """
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
    """
    uses samtools rmdup to quickly eliminate most PCR artifacts.

    writes new $sample.rmd.bam into $results/$sample/$sample.rmdup.bam
    """
    submit = bsub("removedups", P=config['project_id'], verbose=True)

    for sample in config['samples']:
        bam = results / sample / "{sample}.bam".format(sample=sample)
        rmdup = results / sample / "{sample}.rmd.bam".format(sample=sample)
        if rmdup.exists(): continue
        cmd = "samtools rmdup -s {bam} {rmdup}".format(bam=bam.as_posix(),
                                                        rmdup=rmdup.as_posix())
        submit(cmd)


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
    # index all bams >> sample.bam.bai
    jobs = []
    submit = bsub("indexing", P=config['project_id'], verbose=True)
    for bam in results.glob(*.bam):
        bai = Path("{bam}.bai".format(bam=bam.as_posix()))
        if bai.exists(): continue
        cmd = "samtools index {bam}"
        job = submit(cmd)
        jobs.append(job.job_id)



    pass
