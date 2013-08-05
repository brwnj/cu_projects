#!/usr/bin/env python
# encoding: utf-8
"""
Print the track info for given bigwigs.
"""
import os
import argparse

def main(files):
    for file in files:
        file = os.path.basename(file)
        name, ext = os.path.splitext(file)
        parent = "viewPeaks" if "peak" in file else "viewCoverage"
        if "ts" in name.lower():
            stype = "TS"
        elif "nbt" in name.lower():
            stype = "NBT"
        else:
            stype = "UNK"
        view = "PKS" if "peak" in file else "COV"
        inv = "MP" if "MP" in file else "PK"
        if "neg" in file:
            strand = "NEG"
        elif "pos" in file:
            strand = "POS"
        else:
            strand = "B"
        ftype = "bigBed 6" if ext == ".bb" else "bigWig"

        # print ("track {name}\n"
        #         "parent {parent}\n"
        #         "subGroups stype={stype} view={view} inv={inv} strand={strand}\n"
        #         "shortLabel {name}\n"
        #         "longLabel {name}\n"
        #         "bigDataUrl {file}\n"
        #         "color\n"
        #         "type {ftype}\n").format(**locals())

        parent = "mpshifts" if "MP" in file else "pkshifts"
        print ("track {name}\n"
                "parent {parent}\n"
                "subGroups stype={stype} strand={strand}\n"
                "shortLabel {name}\n"
                "longLabel {name}\n"
                "bigDataUrl {file}\n"
                "color 37,52,148\n"
                "type {ftype}\n").format(**locals())

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("files", nargs="+")
    args = p.parse_args()
    main(args.files)
