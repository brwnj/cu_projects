#!/usr/bin/env python
# coding=utf-8
"""
assemble interweaved, paired-end reads; the result of shuffleSequences_fastq.pl.
"""

import logging
import os
import os.path as op
import shutil
import subprocess as sp
import sys
import tempfile as tf
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from itertools import groupby, islice
from toolshed import nopen


# R function wanted a log file name
LOG = None


def verbose_logging(**kwargs):
    logging.info("Record of arguments for %s" % os.path.basename(__file__))
    for k, v in kwargs.items():
        logging.info(">>> %s = %s", k, v)


def readfq(fx):
    with nopen(fx) as fp:
        last = None
        while True:
            if not last:
                for l in fp:
                    if l[0] in '>@':
                        last = l[:-1]
                        break
            if not last: break
            name, seqs, last = last[1:].partition(" ")[0], [], None
            for l in fp:
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+':
                yield name, ''.join(seqs), None
                if not last: break
            else:
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp:
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq):
                        last = None
                        yield name, seq, ''.join(seqs);
                        break
                if last:
                    yield name, seq, None
                    break


def fqreads(fastx):
    total = 0
    for name, seq, qual in readfq(fastx):
        total += 1
    return total


def copyfiles(src, dst, **kwargs):
    cmd = "rsync" + kwargs_to_flag_string(kwargs) + " %s/* %s" % (src, dst)
    return runcmd(cmd)


def copytree(src, dst, symlinks=False, ignore=None):
    """copying to a cifs mount: copy2 will not work."""
    if not op.exists(dst):
        os.makedirs(dst)
        shutil.copystat(src, dst)
    lst = os.listdir(src)
    if ignore:
        excl = ignore(src, lst)
        lst = [x for x in lst if x not in excl]
    for item in lst:
        s = op.join(src, item)
        d = op.join(dst, item)
        if symlinks and op.islink(s):
            if op.lexists(d):
                os.remove(d)
            os.symlink(os.readlink(s), d)
            try:
                st = os.lstat(s)
                mode = stat.S_IMODE(st.st_mode)
                os.lchmod(d, mode)
            except:
                pass # lchmod not available
        elif op.isdir(s):
            copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def runcmd(cmd, log=True):
    """
    >>> import sys
    >>> import logging
    >>> import subprocess as sp
    >>> logging.basicConfig(stream=sys.stdout)
    >>> runcmd("date > /dev/null", False)
    0
    >>> runcmd(["date > /dev/null", "date > /dev/null"], False)
    0
    """
    if isinstance(cmd, basestring):
        if log: logging.debug("Running command: %s", cmd)
        # will raise sp.CalledProcessError on retcode != 0
        sp.check_call(cmd, shell=True)

    else:
        assert isinstance(cmd, list)
        if log: logging.debug("Running commands: %s", cmd)

        # this will be bad if the list of commands is really long
        processes = [sp.Popen(c, shell=True) for c in cmd]
        for i, p in enumerate(processes):
            p.wait()
            if p.returncode != 0:
                raise sp.CalledProcessError(p.returncode, cmd[i])

    return 0


def kwargs_to_flag_string(kwargs, ignore=[]):
    """
    >>> kwargs = {"d":"d_test", "p":"p_test", "sep":","}
    >>> kwargs_to_flag_string(kwargs, ignore=['sep'])
    ' -p p_test -d d_test'
    """
    s = ""
    for k, v in kwargs.items():
        if k in ignore: continue
        s += (" -" if len(k) == 1 else " --") + k + ("" if v is None else (" " + str(v)))
    return s


def run_kmernorm(fastq, **kwargs):
    outfile = fastq.rsplit(".", 1)[0] + ".kmernorm.fastq"
    kwargs['k'] = 19 if not 'k' in kwargs else kwargs['k']
    kwargs['t'] = 80 if not 't' in kwargs else kwargs['t']
    kwargs['c'] = 2 if not 'c' in kwargs else kwargs['c']
    cmd = "kmernorm" + kwargs_to_flag_string(kwargs) + " " + fastq + " > " + outfile
    runcmd(cmd)
    return outfile


def fasta_complexity_filter(fasta, bases='ACTG', threshold=0.01):
    out = fasta.rsplit(".fasta", 1)[0] + ".lowcomp_filter" + ("%.3f" % threshold).lstrip("0") + ".fasta"
    read_count = 0
    kept = 0.

    with open(out, 'wb') as fo:
        for name, seq, qual in readfq(fasta):
            seq_len = float(len(seq))
            read_count += 1
            passed = True

            for b in bases:
                if b not in seq or seq.count(b) / seq_len <= threshold:
                    passed = False

            if passed:
                kept += 1
                fields = []
                print >>fo, "\n".join([">" + name, seq])

    logging.info("low complexity filtering retained %d of %d (%0.3f)", \
                    kept, read_count, kept / read_count * 100)
    return out


def fastq_complexity_filter(fastq, bases='ACTG', threshold=0.01):
    """
    remove low complexity reads (have at least each of the bases specified
    above threshold) and return filtered fastq path.
    """

    def _readfq(fq):
        with nopen(fq) as fh:
            fqclean = (x.strip("\r\n") for x in fh if x.strip())
            while True:
                r = [x for x in islice(fqclean, 8)]
                if not r: raise StopIteration
                assert all(r) and len(r) == 8
                yield r[0], r[1], r[3], r[4], r[5], r[7]

    out = fastq.rsplit(".fastq", 1)[0] + ".lowcomp_filter" + ("%.3f" % threshold).lstrip("0") + ".fasta"
    read_count = 0
    kept = 0.

    with open(out, 'wb') as fo:
        for name1, seq1, qual1, name2, seq2, qual2 in _readfq(fastq):
            # read lengths can vary
            seq1len = float(len(seq1))
            seq2len = float(len(seq2))
            read_count += 2

            passed = True

            for b in bases:
                if b not in seq1 or b not in seq2 \
                        or seq1.count(b) / seq1len <= threshold \
                        or seq2.count(b) / seq2len <= threshold:
                    passed = False

            if passed:
                kept += 2
                fields = [name1, seq1, "+", qual1, name2, seq2, "+", qual2]
                print >>fo, "\n".join(fields)

    logging.info("low complexity filtering retained %d of %d (%0.3f)", \
                    kept, read_count, kept / read_count * 100)
    return out


def spades(**kwargs):
    out = kwargs['o']
    # written for spades 3.0.0; installed in /usr/local/bin
    cmd = "spades.py" + kwargs_to_flag_string(kwargs)
    contigs = op.join(out, "contigs.fasta")
    runcmd(cmd)
    return contigs


def postprocess_spades(fasta, sample, outdir):
    out = op.join(outdir, sample + ".fasta")
    with open(out, 'wb') as fo, nopen(fasta) as fh:
        for line in fh:
            line = line.strip("\r\n")
            if line.startswith(">"):
                line = ">" + sample + "_" + line[1:]
            print >>fo, line
    return out


def assembly_stats(fasta):
    """calculate min_contig_size, max_contig_size, N50, total_bases,
    num_contigs, and gc_percent. writes to log.
    """
    contig_sizes = []
    gc_total = 0

    for name, seq, qual in readfq(fasta):
        contig_sizes.append(len(seq))
        gc_total += seq.count('G') + seq.count('C')

    contigs = len(contig_sizes)

    if contigs > 0:
        contig_sizes.sort(reverse=True)
        maxcontigsize = contig_sizes[0]
        mincontigsize = contig_sizes[-1]
        total_bases = sum(contig_sizes)
        gc_percent = float(gc_total) / total_bases * 100

        if contigs == 1:
            nfifty = contig_sizes[0]
        else:
            testsum = 0
            half = total_bases / 2.
            for size in contig_sizes:
                testsum += size
                if testsum >= half:
                    nfifty = size
                    break

    else:
        total_bases = 0
        mincontigsize = 0
        maxcontigsize = 0
        gc_percent = 0
        nfifty = 0

    logging.info("Stats for FASTA: %s", op.basename(fasta))
    logging.info("Total Bases: %s", total_bases)
    logging.info("Number of Contigs: %s", contigs)
    logging.info("Min Contig Size: %s", mincontigsize)
    logging.info("Max Contig Size: %s", maxcontigsize)
    logging.info("GC Percent: %s", gc_percent)
    logging.info("N50: %s", nfifty)


def read_counts(file, msg):
    count = fqreads(file)
    logging.info("%s: %s", msg, count)


def filter_fasta_by_size(fasta, size=2000):
    path = fasta.rsplit(".", 1)[0]
    passed = "%s.passed.%dbp_min.fasta" % (path, size)
    failed = "%s.failed.%dbp_min.fasta" % (path, size)
    with open(passed, 'w') as passfh, open(failed, 'w') as failfh:
        for name, seq, qual in readfq(fasta):
            if len(seq) > size:
                # may need to word wrap for some applications
                print >>passfh, ">%s\n%s" % (name, seq)
            else:
                print >>failfh, ">%s\n%s" % (name, seq)
    return passed


def send_email(to="scgc@bigelow.org", subject="", message="", attachment=None):
    """ Send a small message via email
        @param to the email address to send the message to
        @param subject the subject line
        @param message a single line message for the email body
        @param attachment a fully qualifed filename of a file to attach (None by default)
        @return the return code for mutt application
        """
    if isinstance(message, list):
        message = " ".join(message)

    muttstub = "echo \"" + message + "\" | mutt"
    if attachment:
        muttstub = muttstub + " -a \"" + attachment + "\""
    cmd = muttstub + " -s \"" + subject + "\" -- " + to
    return runcmd(cmd)


def gzip_all(src, ignore=['.gz']):
    if not '.gz' in ignore: ignore.append('.gz')
    cmds = []

    for f in os.listdir(src):
        f = op.join(src, f)
        if op.isdir(f):
            gzip_all(f, ignore=ignore)
            continue

        append = True
        for extension in ignore:
            if f.endswith(extension):
                append = False
                break

        if append:
            cmds.append("gzip -f %s" % f)

    if len(cmds) > 0:
        runcmd(cmds)


def main(sample, fastq, fasta, output, kmernorm, complexity_filter, email, threads=16):
    verbose_logging(**locals())

    fastq_name = op.basename(fastq)
    fasta_name = op.basename(fasta)

    try:
        tmpdir = tf.mkdtemp("_tmp", "%s_" % sample, tf.tempdir)
        tmpfastq = op.join(tmpdir, fastq_name)
        shutil.copyfile(fastq, tmpfastq)
        tmpfasta = op.join(tmpdir, fasta_name)
        shutil.copyfile(fasta, tmpfasta)

        fq_spades_input = tmpfastq
        fa_spades_input = tmpfasta

        # run kmernorm on the illumina fastq
        if kmernorm:
            kmerfq = run_kmernorm(tmpfastq, k=19, t=80, c=2)
            fq_spades_input = kmerfq

        # complexity filter on fastq and fasta
        if complexity_filter:
            compfilteredfq = lowcomp_filter(fq_spades_input)
            fq_spades_input = compfilteredfq

            compfilteredfa = fasta_complexity_filter(fa_spades_input)
            fa_spades_input = compfilteredfa

            if op.getsize(compfilteredfq) <= 0 or op.getsize(compfilteredfa) <= 0:
                logging.warning("low complexity filtered has removed all of the reads")
                sys.exit(0)

        spades_dir = tmpdir + "/spades"
        # coassembly with spades
        spades_fasta = spades(**{'pacbio':fa_spades_input, 'pe1-12':fq_spades_input,
                                 't':threads, 'o':spades_dir, 'sc':None, 'careful':None})

        # moves spades contigs.fasta and adds sample to header name
        renamed_hdrs = postprocess_spades(spades_fasta, sample, tmpdir)

        read_counts(tmpfastq, "Original Illumina")
        read_counts(tmpfasta, "Original PacBio")
        if kmernorm:
            read_counts(kmerfq, "Digital normalization")
        if complexity_filter:
            read_counts(compfilteredfq, "Complexity filtered Illumina")
            read_counts(compfilteredfa, "Complexity filtered PacBio")

        assembly_stats(renamed_hdrs)
        sizefilteredfa = filter_fasta_by_size(renamed_hdrs, 2000)
        assembly_stats(sizefilteredfa)

    except Exception as e:
        print e.__doc__
        print e.message
        raise

    finally:
        os.remove(tmpfastq)
        os.remove(tmpfasta)
        # gzip all of the files in the temp dir
        gzip_all(tmpdir)
        # copy over the files
        runcmd("cp -R -v {src}/* {dst}".format(src=tmpdir, dst=output))
        # delete the temp working directory
        shutil.rmtree(tmpdir)

        if email:
            send_email(to=email, subject=op.basename(__file__),
                message="finished processing %s; results were copied to %s" % (fastq, output))

        logging.info("Complete.")


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('sample', help="name of sample being processed; used in file naming")
    p.add_argument('fastq', help="interweaved, paired-end reads in fastq format")
    p.add_argument('fasta', help="PACBIO reads to use in co-assembly")
    p.add_argument('output', help="location to store output files")
    p.add_argument('--kmernorm', action='store_true', help="run kmer normalization prior to spades")
    p.add_argument('--complexity-filter', action='store_true', help="filter out low complexity reads before running spades")
    p.add_argument('--email', default="", help="send completion alert")
    p.add_argument('--threads', default=16, type=int, help="threads for spades to utilize")
    args = p.parse_args()

    if not op.exists(args.output):
        os.makedirs(args.output)

    LOG = "%s/%s.log" % (args.output, op.basename(__file__).rstrip(".py"))

    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                        filename=LOG,
                        level=logging.DEBUG,
                        filemode='wb')

    tf.tempdir = tf.gettempdir()
    main(args.sample, args.fastq, args.fasta, args.output, args.kmernorm, \
            args.complexity_filter, args.email, args.threads)