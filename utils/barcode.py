#!/usr/bin/env python
# encoding: utf-8
"""
Analyze potential barcodes of a fastq.
"""
import sys
from collections import Counter

def fastq_reader(fastq):
    

def main(args):
    bc = Counter()
    
    
    for b in reader(args.ta, header="chrom start stop name gene pos count".split()):
        if b['gene'] == ".": continue
        # all ta sites
        cgene.update([b['gene']])
        
        # all that are deemed dire
        if 80 > int(b['pos'].rstrip("%")) > 5:
            cdire.update([b['gene']])
        
        # all that passed cutoff
        if int(b['count']) >= cutoff:
            ccutoff.update([b['gene']])
    return cgene, cdire, ccutoff

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("FASTQ")
    req = p.add_argument_group("required")
    req.add_argument('--length', '-l', type=int, required=True,
                help='length of barcode')
    main(p.parse_args())