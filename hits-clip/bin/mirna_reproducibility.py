#!/usr/bin/env python
# encoding: utf-8
"""
HS5 vs HS27A
HS5 vs hMSC
HS27A vs hMSC
BMEC vs hMSC
BMEC vs HS27A
BMEC vs HS5
BMEC vs HUVEC
MCF7 vs MCF7estr
MCF7 vs BT474
MCF7 vs MDA231
BT474 vs BT474estr
BT474 vs MDA231

python bin/mirna_reproducibility.py -a data/HS5.abundance.noa -b data/HS27A.abundance.noa
python bin/mirna_reproducibility.py -a data/HS5.abundance.noa -b data/hMSC.abundance.noa
python bin/mirna_reproducibility.py -a data/HS27A.abundance.noa -b data/hMSC.abundance.noa
python bin/mirna_reproducibility.py -a data/BMEC.abundance.noa -b data/hMSC.abundance.noa
python bin/mirna_reproducibility.py -a data/BMEC.abundance.noa -b data/HS27A.abundance.noa
python bin/mirna_reproducibility.py -a data/BMEC.abundance.noa -b data/HS5.abundance.noa
python bin/mirna_reproducibility.py -a data/BMEC.abundance.noa -b data/HUVEC.abundance.noa
python bin/mirna_reproducibility.py -a data/MCF7.abundance.noa -b data/MCF7estr.abundance.noa
python bin/mirna_reproducibility.py -a data/MCF7.abundance.noa -b data/BT474.abundance.noa
python bin/mirna_reproducibility.py -a data/MCF7.abundance.noa -b data/MDA231.abundance.noa
python bin/mirna_reproducibility.py -a data/BT474.abundance.noa -b data/BT474estr.abundance.noa
python bin/mirna_reproducibility.py -a data/BT474.abundance.noa -b data/MDA231.abundance.noa
python bin/mirna_reproducibility.py -a data/HELA.abundance.noa -b data/HS27A.abundance.noa
python bin/mirna_reproducibility.py -a data/HELA.abundance.noa -b data/HS5.abundance.noa
python bin/mirna_reproducibility.py -a data/HELA.abundance.noa -b data/hMSC.abundance.noa
"""
import os
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from toolshed import reader
import pandas as pd
import numpy as np

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def get_scaling_factor(abam, bbam, verbose):
    """from two bams, calculate scaling factor for b based on number of mapped reads."""
    if verbose:
        sys.stderr.write(">> Counting reads in A\n")
        #using pysam, count mapped reads in a
        #or just use nopen and samtools view | wc -l
    if verbose:
        sys.stderr.write(">> Counting reads in A\n")
        #using pysam, count mapped reads in b
    sf = float(a) / float(b)
    return sf


def abundance_to_dataframe(a_abundance, b_abundance):
    cases = {}
    for f in [a_abundance, b_abundance]:
        name = os.path.basename(f).split(".")[0]
        cases[name] = {}
        for mirna in reader(f, header="mirname equals abundance".split(), sep=" "):
            # dealing with bad header
            if mirna['mirname'].startswith("hsa"):
                cases[name][mirna['mirname']] = float(mirna['abundance'])
    return pd.DataFrame(cases)


# Polynomial Regression
def polyfit(x, y, degree):
    results = {}

    coeffs = np.polyfit(x, y, degree)

     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = [p(z) for z in x]
    ybar = sum(y)/len(y)
    ssreg = sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = sum([ (yi - ybar)**2 for yi in y])
    results['determination'] = ssreg / sstot

    print results
    return results
    

def main(a, b):
    # imposed on group b
    #sf = get_scaling_factor(args.a, args.b, args.v)
    library_sizes = {'BT474herc':4984545.,'MCF7':17814774.,'MCF7estr':9492187.,
                 'MDA231':14113797.,'BT474':15377810.,'BT474estr':11085221.,
                 'HS27A':17073385.,'HS5':18408107.,'hMSC':10920072.,
                 'BMEC':7885651.,'HUVEC':9805668.,'HELA':17196872.}
    aname = os.path.basename(args.a).split(".")[0]
    bname = os.path.basename(args.b).split(".")[0]
    df = abundance_to_dataframe(args.a, args.b)
    
    # scaling
    df[bname] = df[bname] * (library_sizes[aname] / library_sizes[bname])
    # log2
    df = df.apply(np.log2)
    # correlation coefficient
    c = df[aname].corr(df[bname], method='pearson')
    r2 = c**2
    
    # plot
    fig = Figure(figsize=(6,6))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    # ax.set_title("%s vs %s" % (bname, aname), fontsize=14)
    # ax.set_xlabel(aname, fontsize=12)
    # ax.set_ylabel(bname, fontsize=12)
    ax.set_xlim(left=0, right=25)
    ax.set_ylim(bottom=0, top=25)
    
    # ax.text(24, 1, "R = %.4f" % c, horizontalalignment='right', verticalalignment='bottom')

    # ax.grid(True, linestyle='-', color='0.75')
    ax.scatter(df[aname], df[bname], s=8, color='black')
    
    # Save the generated Scatter Plot to a PNG file.
    canvas.print_figure('%s_vs_%s__rval_%s.png' % (aname, bname, c))


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-a', required=True, help='group a miRNA abundance')
    p.add_argument('-b', required=True, help='group b miRNA abundance')
    args = p.parse_args()
    
    main(args.a, args.b)