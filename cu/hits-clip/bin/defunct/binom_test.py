#!/usr/bin/env python
# encoding: utf-8
"""
"""

import sys
import rpy2
import pandas as pd
import pandas.rpy.common as com
from rpy2.robjects.packages import importr

stats = importr('stats')

def main(args):
    df = pd.io.parsers.read_table(args.highlow, header=0)
    genes = df.columns
    mirs = df.columns
    testing = {}
    for mir in mirs:
        if not mir.startswith("hsa"): continue
        testing[mir] = {}
        # print mir
        for gene in genes:
            if gene.endswith("_notest") or gene.startswith("hsa"): continue
            success = 0
            fail = 0
            for key, row in df.iterrows():
                if row[mir] == 0:
                    # 0 success; 1 fail
                    if row[gene] == 0:
                        success += 1
                    else:
                        fail += 1
                # miR == 0
                else:
                    # 1 success; 0 fail
                    if row[gene] == 1:
                        success += 1
                    else:
                        fail += 1
            # create dataframe for testing
            temp_df = pd.DataFrame({'success': [success], 'fail': [fail]})
            r_matrix = com.convert_to_r_matrix(temp_df)
            p = stats.binom_test(r_matrix)[2][0]
            # print "%s: %f" % (gene, p)
            testing[mir][gene] = p
    testing_df = pd.DataFrame(testing)
    testing_df.to_csv(sys.stdout, sep="\t")

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("highlow", help="series matrix already converted to 1/0")
    main(p.parse_args())