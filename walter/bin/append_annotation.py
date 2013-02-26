#!/usr/bin/env python
"""
Append deseq output onto processed blast output.
"""
import re
import sys
import itertools
from toolshed import reader, header

def make_dict(peaks):
    d = {}
    for line in reader(peaks):
        d[line['trimmedPeak']] = line
    return d

def gtf2dict(gtf):
    header = "chrom source feature start stop score strand frame attributes".split()
    d = {}
    for line in reader(gtf, header=header):
        enst = re.findall(r'transcript_id \"([\w\.]+)\"', line['attributes'])[0].strip()
        ensg = re.findall(r'gene_id \"([\w\.]+)\"', line['attributes'])[0].strip()
        symbol = re.findall(r'gene_name \"([^\"]+)\"', line['attributes'])[0].strip()
        d[enst] = {'ensg':ensg, 'enst':enst, 'symbol':symbol}
    return d

def main(args):
    annotation_dict = make_dict(args.annotated_peaks)
    gtf = gtf2dict(args.gtf)
    blast_header = header(args.blast)
    # only save append these columns onto annotation
    annotation_save = "baseMeanA baseMeanB pval padj sequence".split()
    blast_header.extend(annotation_save)
    # stuff one may want from gtf
    gtf_save = "symbol ensg enst".split()
    blast_header.extend(gtf_save)
    print "\t".join(blast_header)
    for b in reader(args.blast, header=True):
        annotation_lookup = annotation_dict[b['qseqid']]
        # ENST00000576218_-
        enst = b['sseqid'].split("_")[0]
        try:
            gtf_lookup = gtf[enst]
        except KeyError:
            gtf_lookup = {'enst':enst, 'ensg':'-', 'symbol':'-'}
        d = dict(b.items() + annotation_lookup.items() + gtf_lookup.items())
        print "\t".join([d[i] for i in blast_header])

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('blast', help="blast output after filtering for orientation")
    p.add_argument('annotated_peaks', help="deseq output with sequence appended")
    p.add_argument('gtf', help="gtf downloaded straight from ensembl")
    main(p.parse_args())