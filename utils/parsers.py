#! /usr/bin/env python
from itertools import islice, groupby
from toolshed import nopen

def read_fastq(fq):
    with nopen(fq) as fh:
        while True:
            values = list(islice(fh, 4))
            if len(values) == 4:
                id1, seq, id2, qual = values
            elif len(values) == 0:
                raise StopIteration
            else:
                raise EOFError("unexpected end of file")
            assert id1.startswith('@'),\
                    ">> Fastq out of sync at read:\n%s\n" % id1
            assert id2.startswith('+'),\
                    ">> Fastq out of sync at read:\n%s\n" % id1
            assert len(seq) == len(qual),\
                    ">> Sequence and Quality are not the same length \
                    for read:\n%s\n" % id1
            yield id1[1:-1], seq[:-1], qual[:-1]

def read_fasta(fa):
    with nopen(fa) as fh:
        for header, group in groupby(fh, lambda line: line[0] == '>'):
            if header:
                line = group.next()
                name = line[1:].strip()
            else:
                seq = ''.join(line.strip() for line in group)
                yield name, seq

def write_fasta(name, seq, wrap=70):
    print ">%s" % name
    # sequence lines standardized at length of wrap
    print "\n".join([seq[i:i + wrap] for i in range(0, len(seq), wrap)])
