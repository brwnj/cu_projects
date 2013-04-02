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
# stats.binom_test(c(success, fail))
    df = pd.io.parsers.read_table(args.highlow, header=0)
    
    for gene in args.genes.split(","):
        success = 0
        fail = 0
        for key, row in df.iterrows():
            if row[args.miR] == 0:
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
        print "%s: %f" % (gene, p)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("highlow", help="series matrix already converted to 1/0")
    p.add_argument("miR", help="miR to use -- must be same as header")
    p.add_argument("genes", help="comma separated list of genes to test -- must be same as header")
    main(p.parse_args())