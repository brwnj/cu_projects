#!/usr/bin/env python
# coding=utf-8
"""
"""

import yaml
import os.path as op
from bsub import bsub
from glob import glob
from string import Template
from collections import OrderedDict, defaultdict
from argparse import ArgumentParser, RawDescriptionHelpFormatter


def dict_representer(dumper, data):
    return dumper.represent_mapping(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, data.iteritems())


def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))


def load_config(config_file):
    assert op.isfile(config_file), "config file not found"

    yaml.add_representer(OrderedDict, dict_representer)
    yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor)

    with open(config_file) as cf:
        config = yaml.load(cf)
    config['running_jobs'] = defaultdict(list)
    return config


def trim_sequences(config, log):
    """
    trim reads using fastx toolkit.


    config - configuration dictionary
    log - open file handle for verbosity
    """
    trim_cmd = Template(("zcat $fastq | "
                        "fastx_clipper -a $adapter -l $minlength -c -n -v -Q33 | "
                        "fastx_trimmer -Q33 -f 2 | "
                        "gzip -c > $trimmed"))

    # this is added to fastq file name
    pattern = "_trimmed"
    adapter = config['pipeline']['trim']['adapter']
    minlength = config['pipeline']['trim']['minlength']

    submit = bsub("trim", P=config['project_id'], verbose=True)

    for sample, fastq in config['samples'].iteritems():
        assert op.exists(fastq)

        file_path, ext = op.splitext(fastq)
        if ext == ".gz":
            file_path, ext = op.splitext(file_path)

        # trimmed file; single-end specific
        result_file = "%s%s.fastq.gz" % (file_path, pattern)

        if op.exists(result_file):
            # update the fastq path
            config['samples'][sample] = result_file
            print >>log, ">> trimming already complete for", sample
            continue

        # trim
        cmd = trim_cmd.substitute(fastq=fastq,
                                  adapter=adapter,
                                  minlength=minlength,
                                  trimmed=result_file)
        # submit the job
        job = submit(cmd)
        # update config with job ids to wait on
        config['running_jobs'][sample].append(job.job_id)
        # update fastq path in config
        config['samples'][sample] = result_file
        print >>log, ">> trimming started for", sample

   return config


def filter_sequences(config, log):
    """
    Map against rRNA database using bowtie.

    config - configuration dictionary
    log - open file handle for verbosity
    """

    trim_cmd = Template(("zcat $fastq | fastx_clipper -a $adapter -l $minlength "
                        "-c -n -v -Q33 | fastx_trimmer -Q33 -f 2 "
                        "| gzip -c > $trimmed"))

    # this is added to fastq file name
    pattern = "_trimmed"
    adapter = config['pipeline']['trim']['adapter']
    minlength = config['pipeline']['trim']['minlength']

    submit = bsub("trim", P=config['project_id'], verbose=True)

    for sample, fastq in config['samples'].iteritems():
        assert op.exists(fastq)

        file_path, ext = op.splitext(fastq)
        if ext == ".gz":
            file_path, ext = op.splitext(file_path)

        # trimmed file; single-end specific
        result_file = "%s%s.fastq.gz" % (file_path, pattern)
        if op.exists(result_file):
            # update the fastq path
            config['samples'][sample] = result_file
            print >>summary, ">> trimming already complete for", sample
            continue

        # trim
        cmd = trim_cmd.substitute(fastq=fastq,
                                  adapter=adapter,
                                  minlength=minlength,
                                  trimmed=result_file)
        # submit the job
        job = submit(cmd)
        # update config with job ids to wait on
        config['running_jobs'][sample].append(job.job_id)
        # update fastq path in config
        config['samples'][sample] = result_file
        print >>summary, ">> trimming started for", sample

    # print more message to log
    print >>summary, "After Trimming:"
    yaml.dump(config, summary, default_flow_style=False)
    summary.close()
#    return config


# @command('align')
# def align():
#     """
#     align the reads using Novoalign.
#
#     writes $sample.bam to $results/$sample/$sample.bam
#     also writes an alignment summary in the same dir as
#     $sample.alignment_stats.txt
#     """
#     cpus = str(config['align']['cpus'])
#     novoidx = Path(config['align']['index'])
#     assert novoidx.exists(), "index not found"
#
#     submit = bsub("align", P=project_id, n=cpus, R="select[mem>16] rusage[mem=16] span[hosts=1]", verbose=True)
#
#     for sample in config['samples']:
#         fastq = fastqs / '{sample}.fastq.gz'.format(sample=sample)
#
#         # check results dir
#         sample_results = results / sample
#         if not sample_results.exists():
#             sample_results.mkdir()
#
#         bam = results / sample / '{sample}.bam'.format(sample=sample)
#         if bam.exists(): continue
#
#         alignment_summary = results / sample / '{sample}.alignment_stats.txt'.format(sample=sample)
#         align = NOVOALIGN.substitute(index=novoidx.as_posix(),
#                                     fastq=fastq.as_posix(),
#                                     cpus=cpus,
#                                     summary=alignment_summary,
#                                     sample=sample,
#                                     bam=bam)
#         samtools_index = "samtools index {bam}".format(**locals())
#         job = submit(align).then(samtools_index, n=1, R="select[mem>4] rusage[mem=4]")


# @command('bam2bw')
# def bam2bw():
#     """
#     convert all bams to stranded bedgraphs
#     convert those bedgraphs to bigwigs
#
#     from /path/to/something.bam to:
#             /path/to/something_{pos,neg}.bedgraph.gz
#             /path/to/something_{pos,neg}.bw
#     """
#     submit = bsub("bam2bw", P=project_id, verbose=True)
#     symbols, strands = ["+", "-"], ["pos", "neg"]
#     for bam in results.glob('*/*.bam'):
#         base = bam.parent / bam.stem
#
#         p_bedgraph = Path("{base}_pos.bedgraph.gz".format(**locals()))
#         n_bedgraph = Path("{base}_neg.bedgraph.gz".format(**locals()))
#         p_bigwig = Path("{base}_pos.bw".format(**locals()))
#         n_bigwig = Path("{base}_neg.bw".format(**locals()))
#
#         # don't run unnecessarily
#         if p_bedgraph.exists() and n_bedgraph.exists() and \
#                 p_bigwig.exists() and n_bigwig.exists(): continue
#
#         # running the conversions
#         for symbol, strand in zip(symbols, strands):
#             bedgraph = "{base}_{strand}.bedgraph".format(**locals())
#             bigwig = "{base}_{strand}.bw".format(**locals())
#
#             makebg = ("genomeCoverageBed -strand {symbol} -bg -ibam {bam} "
#                         "| bedtools sort -i - > {bedgraph}").format(**locals())
#             makebw = ("bedGraphToBigWig {bedgraph} {chrom_sizes} "
#                         "{bigwig}").format(bedgraph=bedgraph,
#                                             chrom_sizes=chrom_sizes,
#                                             bigwig=bigwig)
#             gzipbg = "gzip -f {bedgraph}".format(**locals())
#
#             job = submit(makebg).then(makebw).then(gzipbg)


def main(config):
    config = load_config(config_file)
    with open(config['results'] + "run_log.out", 'a') as log:
        workflow = config['pipeline'].keys()
        # trim
        if "trim" in workflow:
            config = trim_sequences(config, log)
        # filter
        if "filter" in workflow:
            config = filter_sequences(config, log)
        # map
        # hub
        yaml.dump(config, log, default_flow_style=False)


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('config', help="config file as yaml")
    args = p.parse_args()
    main(args.config)
