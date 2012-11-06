#!/usr/bin/env python
# encoding: utf-8
"""
bedtools intersect output to dexseq_counts.py gff format

bedtools intersect -wb -s -a polya_3utronly.gtf -b refseq.utr3.bed

chr1	hg18_polyaDb	CDS	1158961	1158961	1000.000000	+	0	gene_id "p.126792.1"; transcript_id "p.126792.1"; 	chr1	1158511	1160283	B3GALT6	0	+
chr1	hg18_polyaDb	CDS	1159423	1159423	1000.000000	+	0	gene_id "p.126792.2"; transcript_id "p.126792.2"; 	chr1	1158511	1160283	B3GALT6	0	+
chr1	hg18_polyaDb	CDS	1217268	1217268	1000.000000	+	0	gene_id "p.6339.1"; transcript_id "p.6339.1"; 	chr1	1216853	1217272	SCNN1D	0	+
chr1	hg18_polyaDb	CDS	1236845	1236845	1000.000000	+	0	gene_id "p.54973.1"; transcript_id "p.54973.1"; 	chr1	1236622	1236920	PUSL1	0	+
chr1	hg18_polyaDb	CDS	1236919	1236919	1000.000000	+	0	gene_id "p.126789.1"; transcript_id "p.126789.1"; 	chr1	1236622	1236920	PUSL1	0	+
chr1	hg18_polyaDb	CDS	1366007	1366007	1000.000000	+	0	gene_id "p.64856.1"; transcript_id "p.64856.1"; 	chr1	1364323	1368125	VWA1	0	+
chr1	hg18_polyaDb	CDS	1366007	1366007	1000.000000	+	0	gene_id "p.64856.1"; transcript_id "p.64856.1"; 	chr1	1365030	1368125	VWA1	0	+
chr1	hg18_polyaDb	CDS	1368121	1368121	1000.000000	+	0	gene_id "p.64856.2"; transcript_id "p.64856.2"; 	chr1	1364323	1368125	VWA1	0	+
chr1	hg18_polyaDb	CDS	1368121	1368121	1000.000000	+	0	gene_id "p.64856.2"; transcript_id "p.64856.2"; 	chr1	1365030	1368125	VWA1	0	+

chr1	hg19.gtf.gz	aggregate_gene	11869	14412	.	+	.	gene_id "ENSG00000223972"
chr1	hg19.gtf.gz	exonic_part	11869	11871	.	+	.	transcripts "ENST00000456328"; exonic_part_number "001"; gene_id "ENSG00000223972"
chr1	hg19.gtf.gz	exonic_part	11872	11873	.	+	.	transcripts "ENST00000456328+ENST00000515242"; exonic_part_number "002"; gene_id "ENSG00000223972"
chr1	hg19.gtf.gz	exonic_part	11874	12009	.	+	.	transcripts "ENST00000456328+ENST00000515242+ENST00000518655"; exonic_part_number "003"; gene_id "ENSG00000223972"
chr1	hg19.gtf.gz	exonic_part	12010	12057	.	+	.	transcripts "ENST00000456328+ENST00000515242+ENST00000450305+ENST00000518655"; exonic_part_number "004"; gene_id "ENSG00000223972"
chr1	hg19.gtf.gz	exonic_part	12058	12178	.	+	.	transcripts "ENST00000456328+ENST00000515242+ENST00000518655"; exonic_part_number "005"; gene_id "ENSG00000223972"
chr1	hg19.gtf.gz	exonic_part	12179	12227	.	+	.	transcripts "ENST00000456328+ENST00000515242+ENST00000450305+ENST00000518655"; exonic_part_number "006"; gene_id "ENSG00000223972"
chr1	hg19.gtf.gz	exonic_part	12595	12612	.	+	.	transcripts "ENST00000518655"; exonic_part_number "007"; gene_id "ENSG00000223972"
chr1	hg19.gtf.gz	exonic_part	12613	12697	.	+	.	transcripts "ENST00000518655+ENST00000515242+ENST00000450305+ENST00000456328"; exonic_part_number "008"; gene_id "ENSG00000223972"
chr1	hg19.gtf.gz	exonic_part	12698	12721	.	+	.	transcripts "ENST00000518655+ENST00000515242+ENST00000456328"; exonic_part_number "009"; gene_id "ENSG00000223972"

"%03d" % ( i+1 )
"""
import re
import sys
from toolshed import reader
from collections import defaultdict, OrderedDict

def main(args):
    genes = defaultdict(OrderedDict)
    for i, toks in enumerate(reader(args.intersect, header="chrom s s start stop s strand s attrs s s s genename s s".split())):
        if i == 10:
            break
        polyaname = re.findall(r'transcript_id \"([a-zA-Z-_\/|.0-9]+)\"', toks['attrs'])[0].strip()
        transcript = {}
        transcript[polyaname] = {'chrom':toks['chrom'],'start':toks['start'],'stop':toks['stop'],'strand':toks['strand']}
        genes[toks['genename']].update(transcript)
    # print genes
    
    for gene, transcripts in genes.iteritems():
        # print transcripts
        genestart = genestop = None
        print transcripts
        # get start and stop of the gene
        for polya, toks in transcripts.iteritems():
            if genestart is None:
                genestart = int(toks['start'])
            if genestart > int(toks['start']):
                genestart = int(toks['start'])
            if genestop is None:
                genestop = int(toks['stop'])
            if genestop < int(toks['stop']):
                genestop = int(toks['stop'])
        print genestart
        print genestop

        # print gene line
        # print each polya site
    

    In [1]: import metaseq
    ---------------------------------------------------------------------------
    ImportError                               Traceback (most recent call last)
    <ipython-input-1-3365faa9d779> in <module>()
    ----> 1 import metaseq

    /Users/brownj/devel/metaseq/metaseq/__init__.py in <module>()
          4 from helpers import data_dir, example_filename, nice_colormap, \
          5         gfffeature_to_interval
    ----> 6 from genomic_signal import genomic_signal
          7 import plotutils
          8 import integration

    /Users/brownj/devel/metaseq/metaseq/genomic_signal.py in <module>()
         32 from bx.bbi.bigwig_file import BigWigFile
         33 
    ---> 34 from array_helpers import _array, _array_parallel, _local_coverage, \
         35     _local_coverage_bigwig, _local_count
         36 import filetype_adapters

    /Users/brownj/devel/metaseq/metaseq/array_helpers.py in <module>()
          6 import genomic_signal
          7 import sys
    ----> 8 from rebin import rebin, float_rebin
          9 from helpers import chunker
         10 import filetype_adapters

    ImportError: No module named rebin


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('intersect', help='result file of `bedtools intersect -wb`')
    args = p.parse_args()
    main(args)