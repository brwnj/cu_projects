#!/usr/bin/env python
# encoding: utf-8
"""
"""

from bsub import bsub
from invoke import task
from pathlib import Path
from string import Template
from yaml import dump, load

samples = set('sample1')

project_id = ""
config_path = Path('path/to/yaml')
config = load(open(config_path.as_posix(), 'r'))

NOVOALIGN = Template(("novoalign -d $index -f $fastq -a -o SAM -r A 20 -e 100 -c $cpus -k 2> $summary "
                        "| samtools view -ShuF4 - "
                        "| samtools sort -o - $sample.temp -m 8G "
                        "> $bam"))


@task()
def align():
    """align the reads using Novoalign."""
    cpus = str(config['align']['cpus'])
    submit = bsub("align", P=config['project_id'], n=cpus, R="select[mem>16] rusage[mem=16] span[hosts=1]", verbose=True)
    for sample in config['samples']:
        fastq = data / sample.fastq.gz

        # check results dir
        results = Path(config['results'])
        results = results / sample
        if not results.exists()
            # mkdir

        bam = results / sample / sample.bam
        if bam.exists(): continue
        alignment_summary = results / sample / sample.alignment_stats.txt
        cmd = NOVOALIGN.substitute(index=config['align']['index'],
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
