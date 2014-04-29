#!/usr/bin/env python
# coding=utf-8
"""
parsing maf downloaded from cancer genomics consortium?
"""

import pandas as pd
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def main(data_file, gene, effect, type, find_effect, find_type):
    df = pd.read_table(data_file, compression="gzip" if data_file.ends_with(".gz") else None)
    # groupby gene
    # if

    # track by sample all of the mutations
    # that would give ability to get single sample comutation rates
    # then collapse into


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('data_file', help="maf data file with header")
    data_file_args = p.add_argument_group("required metadata")
    data_file_args.add_argument('--gene', default="Hugo_Symbol", help="column name containing gene ID")
    data_file_args.add_argument('--effect', default="Variant_Class", help="column name containing variant effect")
    data_file_args.add_argument('--type', default="Variant_Type", help="column name containing variant type")
    user_specs_args = p.add_argument_group("parse parameters")
    user_specs_args.add_argument('-e', '--find-effect', default=["Missense_Mutation"], action="append", help="mutation effect to quantify; can be specified multiple times")
    user_specs_args.add_argument('-t', '--find-type', default=["SNP"], action="append", help="mutation type to quantify; can be specified multiple times")
    args = vars(p.parse_args())
    main(**args)
