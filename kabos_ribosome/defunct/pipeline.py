#!/usr/bin/env python
# coding=utf-8
"""
"""

import os
import sys
import yaml
import errno
import logging
import os.path as op
from bsub import bsub
from collections import OrderedDict
from argparse import ArgumentParser, RawDescriptionHelpFormatter


class Sample(object):
    def __init__(self, name, filename, trackcolor="8,69,148", job_id=None,
                    resultpath="results", group=None):
        self.name = name
        self.filename = filename
        self.trackcolor = trackcolor
        self.job_id = job_id
        self.resultpath = resultpath
        self.group = group
        self.files = {'fastq': filename}

    def __str__(self):
        return self.name

    def __repr__(self):
        return "Sample(%s)" % self.name

    @property
    def filepath(self):
        """
        filename minus the file extension
        """
        fp, ext = op.splitext(self.filename)
        if ext == ".gz":
            fp, ext = op.splitext(fp)
        return fp

    @property
    def compressed(self):
        return self.filename.endswith(".gz")

    def newfile(self, k, f):
        """
        Updating filename to the file that needs to be processed next and also
        saving the path in self.files.
        """
        self.files[k] = f
        self.filename = f


def load_config(config_file):
    if not op.isfile(config_file):
        logging.critical("Exiting: configuration file not found.")
        sys.exit(1)

    def _dict_representer(dumper, data):
        return dumper.represent_mapping(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, data.iteritems())

    def _dict_constructor(loader, node):
        return OrderedDict(loader.construct_pairs(node))

    yaml.add_representer(OrderedDict, _dict_representer)
    yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _dict_constructor)

    with open(config_file) as cf:
        config = yaml.load(cf)
    return config


def start_logging(config):
    import time
    logging.basicConfig(filename=config['results'] + "/pipeline.log",
                        format='%(levelname)s:%(message)s', level=logging.DEBUG)
    logging.info("Run: %s", time.strftime("%Y%m%d -- %H:%M:%S"))


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def build_bsub(config, algorithm, **kwargs):
    """
    >>> from bsub import bsub
    >>> config = {'pipeline': {'filter': {1: 'idx', 'bsub': {'P': 'test', 'R': 'span[hosts=1]', 'n': 10}, 'p': 10}}, 'project_id': 'test'}
    >>> b = build_bsub(config, "filter")
    >>> print b.command.replace("logs/", "")
    bsub -e filter.%J.err -J filter -o filter.%J.out -n 10 -P test -R "span[hosts=1]"
    >>> b = build_bsub(config, "filter", **{'w':10010})
    >>> print b.command.replace("logs/", "")
    bsub -e filter.%J.err -J filter -o filter.%J.out -n 10 -P test -R "span[hosts=1]" -w "done(10010)"
    """
    try:
        pid = config['project_id']
    except KeyError:
        # this is required
        logging.critical("Define a Project ID (project_id) in the config")
        sys.exit(1)

    try:
        # args as defined in config:pipeline:algorithm:bsub
        config_kwargs = config['pipeline'][algorithm]['bsub']
        # overwrite existing with new
        config_kwargs.update(kwargs)
    except KeyError:
        # LSF reservations not defined in config
        if not kwargs:
            return bsub(algorithm, P=pid, verbose=True)
        config_kwargs = kwargs

    # fix wait syntax
    if 'w' in config_kwargs.keys():
        config_kwargs['w'] = '"exit({i},0)"'.format(i=config_kwargs['w'])
    if not 'P' in config_kwargs.keys():
        config_kwargs['P'] = pid

    # args to strings
    for k, v in config_kwargs.items():
        if isinstance(v, int):
            config_kwargs[k] = str(v)

    return bsub(algorithm, verbose=True, **config_kwargs)


def options_to_string(d):
    """
    >>> d = {1:'bowtie_idx', 'm': 1, 'quiet': True}
    >>> options_to_string(d)
    ' -m 1 --quiet bowtie_idx'
    """
    s = ""
    p = []
    for k, v in d.items():
        # handle positional args
        if isinstance(k, int):
            p.insert(k, v)
        elif k != "bsub":
            s += (" --" if len(k) > 1 else " -") + k + \
                        ("" if v is True else (" " + str(v)))
    s += " " + " ".join(p)
    return s


def submit(cmd, sample, config, algorithm, result_file):
    """
    Check if result_file exists or submit job. Updates sample.filename and
    sample.files[algorithm] = result_file

    cmd - string to be executed
    sample - sample object
    config - config dictionary
    algorithm - the current step in the pipeline
    result_file - file being generated by algorithm
    """
    if op.exists(result_file):
        sample.newfile(algorithm, result_file)
        logging.info("%s complete for %s", algorithm, sample)
    else:
        kwargs = {'w': sample.job_id} if sample.job_id else {}
        submit = build_bsub(config, algorithm, **kwargs)
        job = submit(cmd)
        sample.job_id = int(job)
        sample.newfile(algorithm, result_file)
        logging.info("%s queued for %s <%s>", algorithm, sample, sample.job_id)


def trim_sequences(samples, config, pattern="_trimmed"):
    """
    Trim reads using fastx toolkit. Updates current fastq path for each sample.

    samples - list of sample objects
    config - configuration dictionary
    pattern - string appended to resultant fastq
    """
    algorithm = "trim"
    adapter = config['pipeline'][algorithm]['adapter']
    minlength = config['pipeline'][algorithm]['minlength']

    for sample in samples:
        result_file = "%s%s.fastq.gz" % (sample.filepath, pattern)
        stats_file = "%s/%s_stats.txt" % (sample.resultpath, algorithm)
        sample.files['%s_stats' % algorithm] = stats_file

        cmd = ("{cat} {fastq} "
                "| fastx_clipper -a {adapter} -l {minlength} -c -n -v -Q33 2>> {stats}"
                "| fastx_trimmer -Q33 -f 2 2>> {stats}"
                "| gzip -c > {result}").format(cat="zcat" if sample.compressed else "cat",
                                                fastq=sample.filename,
                                                adapter=adapter,
                                                minlength=minlength,
                                                stats=stats_file,
                                                result=result_file)

        submit(cmd, sample, config, algorithm, result_file)


def filter_sequences(samples, config, pattern="_filtered"):
    """
    Map against rRNA database using bowtie. Updates current fastq path for each
    sample.

    samples - list of sample objects
    config - configuration dictionary
    pattern - string appended to resultant fastq
    """
    algorithm = "filter"
    opts = options_to_string(config['pipeline'][algorithm])

    for sample in samples:
        result_file = "%s%s.fastq.gz" % (sample.filepath, pattern)
        temp_file = "%s_temp.fastq" % (sample.filepath)
        stats_file = "%s/%s_stats.txt" % (sample.resultpath, algorithm)
        sample.files['%s_stats' % algorithm] = stats_file

        filt = ("{cat} {fastq} "
                "| bowtie --un {result} {opts} - "
                "> /dev/null "
                "2> {stats}").format(cat="zcat" if sample.compressed else "cat",
                                     fastq=sample.filename,
                                     result=temp_file,
                                     opts=opts,
                                     stats=stats_file)
        gzip = "gzip {fastq}".format(fastq=temp_file)
        rename = "mv {temp}.gz {result}".format(temp=temp_file,
                                                result=result_file)

        submit(filt, sample, config, algorithm, result_file)
        submit(gzip, sample, config, "gzip", result_file)
        submit(rename, sample, config, "rename", result_file)


def tophat(samples, config, pattern="_tophat"):
    """
    Align reads using tophat.

    samples - list of sample objects
    config - configuration dictionary
    pattern - string appended to resultant bam
    """
    algorithm = "tophat"
    opts = options_to_string(config['pipeline'][algorithm])
    for sample in samples:
        result_file = "%s/%s%s.bam" % (sample.resultpath, sample, pattern)
        temp_file = "%s/accepted_hits.bam" % sample.resultpath
        stats_file = "%s/align_summary.txt" % sample.resultpath
        sample.files['%s_stats' % algorithm] = stats_file

        align = ("tophat {opts} --output-dir {out} "
                    "{fastq}").format(opts=opts,
                                      out=sample.resultpath,
                                      fastq=sample.filename)
        rename = "mv {temp_bam} {sample_bam}".format(temp_bam=temp_file,
                                                     sample_bam=result_file)
        index = "samtools index {bam}".format(bam=result_file)

        # should return bsub object and use .then()
        submit(align, sample, config, algorithm, result_file)
        submit(rename, sample, config, "rename", result_file)
        submit(index, sample, config, "index", result_file)


# def bam2bw():
#     """
#     """
#     for sample in samples:
#         for bam in files:
#
#
#     symbols, strands = ["+", "-"], ["pos", "neg"]
#     for bam in results.glob('*/*.bam'):
#         base = bam.parent / bam.stem
#
#         p_bedgraph = Path("{base}_pos.bedgraph.gz".format(**locals()))
#         n_bedgraph = Path("{base}_neg.bedgraph.gz".format(**locals()))
#         p_bigwig = Path("{base}_pos.bw".format(**locals()))
#         n_bigwig = Path("{base}_neg.bw".format(**locals()))
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


def get_samples(config):
    samples = []

    for name, kwargs in config['samples'].items():
        # all of the initial fastqs should exist
        try:
            filename = kwargs.pop('fastq')
        except KeyError:
            logging.critical("Exiting. 'fastq' must be defined for %s", name)
            sys.exit(1)
        if not op.exists(filename):
            logging.critical("Path for %s isn't valid (%s).", name, filename)
            sys.exit(1)

        resultpath = "%s/%s" % (config['results'], name)
        # should also be able to write here
        try:
            mkdir_p(resultpath)
        except OSError:
            logging.critical("Exiting. Can't create directory %s.", resultpath)
            sys.exit(1)
        kwargs['resultpath'] = resultpath

        samples.append(Sample(name, filename, **kwargs))

    return samples


def main(config):
    config = load_config(config)
    start_logging(config)

    workflow = config['pipeline'].keys()
    samples = get_samples(config)

    # could be generalized where config key is method name

    if "trim" in workflow:
        trim_sequences(samples, config)

    if "filter" in workflow:
        filter_sequences(samples, config)

    if "tophat" in workflow:
        tophat(samples, config)

    # hub
    # yaml.dump(config, default_flow_style=False)


if __name__ == '__main__':
    import doctest
    if doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE).failed == 0:
        p = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
        p.add_argument('config', help="Configuration file as yaml.")
        args = p.parse_args()
        main(args.config)