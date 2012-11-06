#!/usr/bin/env python
# encoding: utf-8
"""
HS5 vs HS27A
HS5 vs HS27A vs hMSC
HS5 vs hMSC
HS27A vs hMSC
BMEC vs hMSC
BMEC vs HS27A vs HS5
MCF7 vs MCF7estr
BT474 vs BT474estr
MCF7 vs BT474
MCF7 vs MDA231
BT474 vs MDA231
"""
import argparse
import os
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from toolshed import reader
from itertools import izip
import pandas as pd
import numpy as np

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def get_args():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--abundance', action='append', required=True, help='miRNA abundance')
    p.add_argument('--microarray', required=True, help='hs27a vs hs5 microarray data')
    args = p.parse_args()
    return args


def foldchange(a, b):
    amean = np.mean(a)
    bmean = np.mean(b)
    fc = 0.
    if amean != 0. and bmean != 0.:
        if amean > bmean:
            # just taking absolute value, so no * -1
            fc = np.power(2, np.log2(amean) - np.log2(bmean))
        if amean < bmean:
            fc = np.power(2, np.log2(bmean) - np.log2(amean))
    else:
        if amean > bmean:
            fc = amean
        if amean < bmean:
            fc = bmean
    return fc


def junk_to_dataframe(abundance_files, microarray):
    cases = {}
    abundance = {}
    for f in abundance_files:
        name = os.path.basename(f).split(".")[0]
        abundance[name] = {}
        for mirna in reader(f, header="mirname equals abundance".split(), sep=" "):
            # inefficient. dealing with stupid header
            if mirna['mirname'].startswith("hsa"):
                abundance[name][mirna['mirname']] = float(mirna['abundance'])
                
    cases['hitsclip'] = {}
    cases['microarray'] = {}
    for mirna in reader(microarray, header="name qvalue hs27a hs27aerror hs5 hs5error direction foldchange".split(), sep="\t"):
        # deal with yet another weird header
        if mirna['name'].startswith("hsa"):
            cases['microarray'][mirna['name']] = float(mirna['foldchange'])
    mirnas = cases['microarray'].keys()
    for mirna in mirnas:
        try:
            cases['hitsclip'][mirna] = foldchange(abundance['HS27A'].get(mirna), abundance['HS5'].get(mirna))
        except:
            mirna = mirna.rstrip("*")
            try:
                cases['hitsclip'][mirna] = foldchange(abundance['HS27A'].get(mirna + "-5p"), abundance['HS5'].get(mirna + "-5p"))
            except:
                try:
                    cases['hitsclip'][mirna] = foldchange(abundance['HS27A'].get(mirna + "-3p"), abundance['HS5'].get(mirna + "-3p"))
                except:
                    continue
    return pd.DataFrame(cases)


def main():
    args = get_args()
    df = junk_to_dataframe(args.abundance, args.microarray)
    
    # correlation coefficient
    #c = df[aname].corr(df[bname], method='pearson')
    # r2 = c**2
    
    df = df.apply(np.log2)
    
    # plot
    fig = Figure(figsize=(10,10))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    # ax.set_title("HS5: Microarray vs Hits-Clip", fontsize=14)
    ax.set_xlabel("Hits-Clip", fontsize=12)
    ax.set_ylabel("Microarray", fontsize=12)
    ax.scatter(df['hitsclip'], df['microarray'], s=8, color='black')
    c = df['hitsclip'].corr(df['microarray'], method='pearson')
    ax.text(0.95, 0.1, "R = %.4f" % c, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    
    # ax = fig.add_subplot(422)
    #     ax.set_title("HS27A: Microarray vs Hits-Clip", fontsize=14)
    #     ax.set_xlabel("Hits-Clip", fontsize=12)
    #     ax.set_ylabel("Microarray", fontsize=12)
    #     ax.scatter(df['HS27A'], df['HS27A_m'], s=20, color='tomato')
    #     
    #     c = df['HS27A'].corr(df['HS27A_m'], method='pearson')
    #     ax.text(0.95, 0.1, "R = %.4f" % c, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    #     
    #     df = df.apply(np.log2)
    #     
    #     ax = fig.add_subplot(423)
    #     ax.set_title("HS5: Microarray vs Hits-Clip (log2)", fontsize=14)
    #     ax.set_xlabel("Hits-Clip", fontsize=12)
    #     ax.set_ylabel("Microarray", fontsize=12)
    #     ax.scatter(df['HS5'], df['HS5_m'], s=20, color='tomato')
    #     
    #     c = df['HS5'].corr(df['HS5_m'], method='pearson')
    #     ax.text(0.95, 0.1, "R = %.4f" % c, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    #     
    #     ax = fig.add_subplot(424)
    #     ax.set_title("HS27A: Microarray vs Hits-Clip (log2)", fontsize=14)
    #     ax.set_xlabel("Hits-Clip", fontsize=12)
    #     ax.set_ylabel("Microarray", fontsize=12)
    #     ax.scatter(df['HS27A'], df['HS27A_m'], s=20, color='tomato')
    #     
    #     c = df['HS27A'].corr(df['HS27A_m'], method='pearson')
    #     ax.text(0.95, 0.1, "R = %.4f" % c, horizontalalignment='right', verticalalignment='center', transform=ax   .transAxes)
    #     
    #     fig.subplots_adjust(hspace=0.4)
    #     
    # Save the generated Scatter Plot to a PNG file.
    canvas.print_figure('micro_correlation.png', dpi=500)


if __name__ == "__main__":
    main()