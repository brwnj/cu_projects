#!/usr/bin/env python
# encoding: utf-8
"""
Group AA sequences into families with the ability to collapse sequences based
on a specified number of sequence mismatches.

Usage:
    groupby_aa IMGT_SEQUENCES [--mismatches=N]
    groupby_aa (-h | --help)
    groupby_aa --version

Examples:
    groupby_aa 5_AA-sequences_4_TGGTCA_5_150813.txt
    groupby_aa 5_AA-sequences_4_TGGTCA_5_150813.txt --mismatches 2

Options:
    -h --help           Show this screen.
    --version           Show version.
    --mismatches NUM    Number of allowed mismatches when grouping CDR3
                        sequences within a family [default: 0].
"""
import pandas as pd
from docopt import docopt

def compression_level(file_name):
    return "gzip" if file_name.endswith("gz") else None

def main(aa_sequences, mismatches):
    # import into a table
    compression = compression_level(aa_sequences)
    df = pd.read_table(aa_sequences, index_col=0, compression=compression)
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
    header = ['V-GENE','J-GENE','D-GENE','CDR3']
    print "\t".join(header)
    # no further collapsing of the CDR3 sequences
    if mismatches == 0:
        # adding CDR3 to groupby collapses them into uniques only
        for (v, j, d, cdr3), (dframe) in df.groupby(['V-GENE and allele','J-GENE and allele','D-GENE and allele','CDR3-IMGT']):
            print "\t".join([v,j,d,cdr3])
    else:
        for (v, j, d), (dframe) in df.groupby(['V-GENE and allele','J-GENE and allele','D-GENE and allele']):
            # remove duplicate cdr3 sequences
            dframe.drop_duplicates('CDR3-IMGT', inplace=True)
            # if len(df) == 1:
                # print ...
        """
        Homsap IGHV1-2*02 F     Homsap IGHJ4*02 F       Homsap IGHD7-27*01 F    ARALPGDAMGGLHY
        ARALPGDAMGGLHY
        ARALTEDAMGGLHY
        ARALTGAAMGGLHY
        ARALTGDAMGGLHY
        ARALTGDAMGGLRY
        ARALTGDAMSGLHY
        ARALTGGAMGGLHY
        ARALTGNAMGGLHY
        
        Homsap IGHV1-2*02 F     Homsap IGHJ5*02 F       Homsap IGHD4-11*01 ORF  AEETRFTITTSFDP
        AEETRFTITTSFDP
        AGETRFTVTTSFDP
        AGGTRFTVTTSFDP
        AIETMFTVTTSVDP
        AKEPRFTVPTSFDP
        AKEPRFTVTTSFDP
        AKETRFTVTTSFDP
        ARATRVTVTTPFDP
        ARDIRFTFTTSIDP
        ARDTRFTVTTSVDP
        AREARFTVTTSFDP
        AREPRFTVTTSFDP
        ARESRFTVTTSFDP
        ARETGFTVTTSFDP
        ARETKFTVTTSFDP
        ARETQFTVTTSFHP
        ARETRCTVTTSFDP
        ARETRFSVTTSFDP
        ARETRFTVATSFDP
        ARETRFTVPTSFDP
        ARETRFTVSTSFDP
        ARETRFTVSTSLDA
        ARETRFTVTNSFDP
        ARETRFTVTSSFDA
        ARETRFTVTSSFDP
        ARETRFTVTTPFDP
        ARETRFTVTTSCDL
        ARETRFTVTTSFDL
        ARETRFTVTTSFDP
        ARETRFTVTTSFGP
        ARETRFTVTTSIDP
        ARETRLTVTTSFDP
        ARETRSTVTTSFDP
        ARETRVTVTTSFDP
        ARGTRFTVTTSFDP
        ARQTRFTVTTSFDP
        ERETRFTVTTSFDP
        SRETRFTVTTSFDP
        TRETRFTVTTSFDP
        """

if __name__ == '__main__':
    args = docopt(__doc__, version='GroupBy AA 0.0.1')
    # docopt is garbage...
    main(args['IMGT_SEQUENCES'], int(args['--mismatches']))
