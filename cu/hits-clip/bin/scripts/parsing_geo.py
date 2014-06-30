#!/usr/bin/env python
# encoding: utf-8
"""
group miRs into high low from transposed series-matrix.
"""
import sys
import pandas as pd

def main(args):
    df = pd.io.parsers.read_csv(args.series_matrix, sep="\t")
    ilmn_ids = df.columns
    for ilmn_id in ilmn_ids:
        if ilmn_id.endswith("_notest"): continue
        df = df.sort(ilmn_id)
        for i, (k, v) in enumerate(df[ilmn_id].iteritems()):
            if i < (len(df) / 2):
                df[ilmn_id][k] = 0
            else:
                df[ilmn_id][k] = 1
    df.to_csv(sys.stdout, sep="\t", index=False)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("series_matrix", help="transposed series matrix")
    main(p.parse_args())