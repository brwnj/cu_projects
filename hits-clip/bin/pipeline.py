#!/usr/bin/env python
# encoding: utf-8
"""
"""

from pathlib import Path

from bsub import bsub
from glob import glob
from string import Template
from yaml import dump, load
from commandr import command, Run
from os.path import basename, exists, isdir, isfile, splitext

config_path = '/vol1/home/brownj/projects/hits-clip/bin/config.yaml'
assert isfile(config_path), "config file not found"
config = load(open(config_path))

fastqs = Path(config['fastqs'])
assert fastqs.is_dir()
results = config['results']
assert exists(results)
chrom_sizes = config['chrom_sizes']
assert isfile(chrom_sizes)
project_id = config['project_id']

TRIM = Template(("python {script} --reverse-complement -a {adapter} "
                    "$untrimmed | gzip -c > $trimmed").format(\
                        script=config['trim_script'],
                        adapter=config['adapter_sequence']))

NOVOALIGN = Template(("novoalign -d $index -f $fastq -a -o SAM -r A 20 -e 100 -c $cpus -k 2> $summary "
                        "| samtools view -ShuF4 - "
                        "| samtools sort -o - $sample.temp -m 8G "
                        "> $bam"))

@command('trimsequence')
def trim_sequences():
    """
    trim reads by design using script that checks for sequence.

    writes new fastq will same name as previous
    renames untrimmed fastq $sample.untrimmed.fastq.gz
    """
    submit = bsub("trim", P=project_id, verbose=True)
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


@command('align')
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

    submit = bsub("align", P=project_id, n=cpus, R="select[mem>16] rusage[mem=16] span[hosts=1]", verbose=True)

    for sample in config['samples']:
        fastq = fastqs / '{sample}.fastq.gz'.format(sample=sample)

        # check results dir
        sample_results = results / sample
        if not sample_results.exists():
            sample_results.mkdir()

        bam = results / sample / '{sample}.bam'.format(sample=sample)
        if bam.exists(): continue

        alignment_summary = results / sample / '{sample}.alignment_stats.txt'.format(sample=sample)
        align = NOVOALIGN.substitute(index=novoidx.as_posix(),
                                    fastq=fastq.as_posix(),
                                    cpus=cpus,
                                    summary=alignment_summary,
                                    sample=sample,
                                    bam=bam)
        samtools_index = "samtools index {bam}".format(**locals())
        job = submit(align).then(samtools_index, n=1, R="select[mem>4] rusage[mem=4]")


# @task()
# def indexbams():
#     """
#     create and index for any bam in the results directory.
#     """
#     submit = bsub("index", P=config['project_id'], verbose=True)
#     for bam in results.glob('*/*.bam'):
#         bai = Path("{bam}.bai".format(bam=bam.as_posix()))
#         if bai.exists(): continue
#         cmd = "samtools index {bam}".format(bam=bam.as_posix())
#         job = submit(cmd)


@command('removedups')
def removedups():
    """
    uses samtools rmdup to quickly eliminate most PCR artifacts.

    writes new $sample.rmd.bam into $results/$sample/$sample.rmdup.bam
    """
    submit = bsub("removedups", P=project_id, verbose=True)

    for sample in config['samples']:
        bam = results / sample / "{sample}.bam".format(sample=sample)
        rmdup = results / sample / "{sample}.rmd.bam".format(sample=sample)
        if rmdup.exists(): continue
        samtools_rmdup = "samtools rmdup -s {bam} {rmdup}".format(
                                            bam=bam.as_posix(),
                                            rmdup=rmdup.as_posix())
        samtools_index = "samtools index {rmdup}".format(rmdup=rmdup.as_posix())

        job = submit(samtools_rmdup).then(samtools_index)


@command('bam2bw')
def bam2bw():
    """
    convert all bams to stranded bedgraphs
    convert those bedgraphs to bigwigs

    from /path/to/something.bam to:
            /path/to/something_{pos,neg}.bedgraph.gz
            /path/to/something_{pos,neg}.bw
    """
    submit = bsub("bam2bw", P=project_id, verbose=True)
    symbols, strands = ["+", "-"], ["pos", "neg"]
    for bam in results.glob('*/*.bam'):
        base = bam.parent / bam.stem

        p_bedgraph = Path("{base}_pos.bedgraph.gz".format(**locals()))
        n_bedgraph = Path("{base}_neg.bedgraph.gz".format(**locals()))
        p_bigwig = Path("{base}_pos.bw".format(**locals()))
        n_bigwig = Path("{base}_neg.bw".format(**locals()))

        # don't run unnecessarily
        if p_bedgraph.exists() and n_bedgraph.exists() and \
                p_bigwig.exists() and n_bigwig.exists(): continue

        # running the conversions
        for symbol, strand in zip(symbols, strands):
            bedgraph = "{base}_{strand}.bedgraph".format(**locals())
            bigwig = "{base}_{strand}.bw".format(**locals())

            makebg = ("genomeCoverageBed -strand {symbol} -bg -ibam {bam} "
                        "| bedtools sort -i - > {bedgraph}").format(**locals())
            makebw = ("bedGraphToBigWig {bedgraph} {chrom_sizes} "
                        "{bigwig}").format(bedgraph=bedgraph,
                                            chrom_sizes=chrom_sizes,
                                            bigwig=bigwig)
            gzipbg = "gzip -f {bedgraph}".format(**locals())

            job = submit(makebg).then(makebw).then(gzipbg)


@command('genomedata')
def genomedata():
    """
    creates genomedata archive. if the archive already exists, adds tracks
    until all project bedgraphs have been added.
    """
    from genomedata import Genome

    genomedata_path = config['genomedata']

    if isdir(genomedata_path):
        # get current track list
        genome = Genome(genomedata_path)
        existing_tracks = genome.tracknames_continuous
        genome.close()

        # iterate over existing bedgraphs
        to_add = {}
        for bg in glob(results + "/*/*.bedgraph.gz"):
            trackname, ext = basename(bg).rsplit(".bedgraph.gz", 1)
            if trackname in existing_tracks: continue
            to_add[trackname] = bg

        assert len(to_add) > 0, "no new data to load"

        # open the genomedata archive
        submit_open = bsub("open_data", P=project_id)
        openarchive = "genomedata-open-data {genomedata_path} {tracks}".format(
                            genomedata_path=genomedata_path,
                            tracks=" ".join(to_add.keys()))
        job = submit_open(openarchive)

        # each job is dependent on the previous to finish
        # load each of the tracks into the archive
        for track, bedgraph in to_add.iteritems():
            submit_load = bsub("load_data", P=project_id, w=str(job.job_id))
            loaddata = "zcat {bedgraph} | genomedata-load-data {genomedata_path} {track}".format(**locals())
            job = submit_load(loaddata)

        submit_close = bsub("close_data", P=project_id, w=str(job.job_id))
        closearchive = "genomedata-close-data {genomedata_path}".format(**locals())
        submit_close(closearchive)

    else:
        # create entire archive
        pass


if __name__ == '__main__':
    Run()
