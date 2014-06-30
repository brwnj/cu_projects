#!/usr/bin/env python
# encoding: utf-8
"""
"""
import os
import re
import sys
import argparse

KEY = {"PK86-7":"PK63","PK93-1":"PK64","PK93-3":"PK65","PK100-3":"PK66",
        "PKm+c":"PK67","PKm+e2":"PK68","PKnbt29":"PK69","PKnbt39":"PK70",
        "PKnbt89":"PK71","PKts21":"PK72","PKts28":"PK73","PKts57":"PK74",
        "PK86-7a":"PK75","PK89-5a":"PK76","PK89-9a":"PK77","PK93-1a":"PK78",
        "PK93-3a":"PK79","PK100-3a":"PK80","PKNbt89":"PK81","PKNbt94":"PK82",
        "PKNbt102":"PK83","PKTs57":"PK84","PKTs61":"PK85","PKTs68":"PK86"}

def main(files):
    for f in files:
        key = re.split('\.|_', f)[0]
        try:
            newname = f.replace(key, KEY[key])
            os.rename(f, newname)
        except KeyError:
            print >>sys.stderr, ">> skipped", f
            continue

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("files", nargs="+")
    args = p.parse_args()
    main(args.files)
