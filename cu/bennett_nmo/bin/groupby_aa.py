#!/usr/bin/env python
# encoding: utf-8
"""
Group AA sequences into OTUs, utilizing cd-hit. cd-hit must be in PATH.
"""

#"Clustering of highly homologous sequences to reduce thesize of large protein database", Weizhong Li, Lukasz Jaroszewski & Adam Godzik. Bioinformatics, (2001) 17:282-283
#"Tolerating some redundancy significantly speeds up clustering of large protein databases", Weizhong Li, Lukasz Jaroszewski & Adam Godzik. Bioinformatics, (2002) 18:77-82

import os
import sys
import argparse
import pandas as pd
import subprocess as sp
from toolshed import nopen
from itertools import groupby

def compression_level(file_name):
    return "gzip" if file_name.endswith("gz") else None

def readfa(fa):
    with nopen(fa) as fh:
        for header, group in groupby(fh, lambda line: line[0] == '>'):
            if header:
                line = group.next()
                name = line[1:].strip()
            else:
                seq = ''.join(line.strip() for line in group)
                yield name, seq

def cluster_proteins(sequences, threshold):
    """
    sequences   pd.Series of AA sequences

    returns list of sequences.
    """
    # create fasta of sequences
    input_fasta = open("cdhit_in.fasta", "wb")
    output_fasta = "cdhit_out.fasta"
    # generate fasta of sequences
    for i, seq in enumerate(sequences):
        print >>input_fasta, ">seq_{id}\n{sequence}".format(id=i, sequence=seq)
    input_fasta.close()
    # run fasta through cd-hit
    cmd = "cd-hit -i {input} -o {output} -c {threshold}".format(input=input_fasta.name,
                                                                output=output_fasta,
                                                                threshold=threshold)

    p = sp.Popen(cmd, stderr=sp.PIPE, stdout=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    # parse output for remaining sequences
    clustered_proteins = [seq for name, seq in readfa(output_fasta)]
    # cleanup files from cd-hit call
    os.remove(input_fasta.name)
    os.remove(output_fasta)
    os.remove("{output}.clstr".format(output=output_fasta))
    # return list of remaining sequences
    return clustered_proteins

def main(aa_sequences, threshold, column_name):
    # import into a table
    compression = compression_level(aa_sequences)
    df = pd.read_table(aa_sequences, index_col=0, compression=compression)
    assert column_name in df.columns
    # df = pd.read_table("5_AA-sequences_4_TGGTCA_5_150813.txt", index_col=0)
    # drop unproductive translations
    df = df[df.Functionality == "productive"]
    # only take one suggestion from the imgt translation for each v j and d
    # some have multiple possibilities separated by commas
    vgene = df['V-GENE and allele'].apply(lambda x: pd.Series(x.split(",")[0]))
    jgene = df['J-GENE and allele'].apply(lambda x: pd.Series(x.split(",")[0]))
    # dtype appears as float in testing
    dgene = df['D-GENE and allele'].astype("string").apply(lambda x: pd.Series(x.split(",")[0]))
    # remove existing entries
    del df['V-GENE and allele']
    del df['J-GENE and allele']
    del df['D-GENE and allele']
    # add everything back in
    df = df.join(vgene)
    df = df.join(jgene)
    df = df.join(dgene)
    header = ['V-GENE', 'J-GENE', 'D-GENE', column_name]
    print "\t".join(header)
    # no further collapsing of the CDR3 sequences
    if threshold == 1.0:
        # adding sequence to groupby collapses them into uniques only
        for (v, j, d, sequence), (dframe) in df.groupby(['V-GENE and allele', 'J-GENE and allele', 'D-GENE and allele', column_name]):
            print "\t".join([v, j, d, sequence])
    else:
        for (v, j, d), (dframe) in df.groupby(['V-GENE and allele', 'J-GENE and allele', 'D-GENE and allele']):
            # remove duplicate sequences
            dframe.drop_duplicates(column_name, inplace=True)
            if len(dframe) == 1:
                try:
                    print "\t".join([v, j, d, dframe[column_name].values[0]])
                except TypeError:
                    # no sequence present; "NaN" in dataframe; ignore
                    continue
            else:
                # cluster the sequences based on threshold
                sequences = cluster_proteins(dframe[column_name], threshold)
                for seq in sequences:
                    print "\t".join([v, d, j, seq])

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('IMGT_SEQUENCES', help="Sequences obtained from HighV-QUEST")
    p.add_argument('-c', '--column', default="CDR3-IMGT",
            help="Column name of which to base sequence clustering")
    p.add_argument('-t', '--identity-threshold', type=float, default=0.90,
            help="Sequence Identity threshold: the number of identical amino \
            acids in alignment divided by the full length of the shorter \
            sequence.")
    args = p.parse_args()
    if 1 > args.identity_threshold < 0.65:
        print >>sys.stderr, "Identity threshold can be between 1 and 0.65, inclusive"
        sys.exit(1)
    main(args.IMGT_SEQUENCES, args.identity_threshold, args.column)
