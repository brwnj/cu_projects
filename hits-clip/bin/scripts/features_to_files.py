#!/usr/bin/env python
# encoding: utf-8
"""
Converts Cytoscape edge features to separate feature files for Cytoscape
import.

Takes:

hsa-miR-34a-5p (7) TMSB4X = cds
hsa-miR-1290 (7) TMSB4X = cds
hsa-miR-2682-5p (7) TMSB4X = cds

And makes:

featureCDS
hsa-miR-34a-5p (7) TMSB4X = True
hsa-miR-1290 (7) TMSB4X = True
hsa-miR-2682-5p (7) TMSB4X = True

featureUTR5
hsa-miR-34a-5p (7) TMSB4X = True
hsa-miR-1290 (7) TMSB4X = True
hsa-miR-2682-5p (7) TMSB4X = True

etc.
"""

import argparse
import os
import sys
from toolshed import nopen

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"

def main():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('attrs', help='Cytoscape attributes file of overlapping features')
    args = p.parse_args()
    
    attrsfile = nopen(args.attrs)
    # remove header
    attrsfile.readline()
    
    uniqid = os.path.basename(args.attrs).split(".", 1)[0]
    
    previous_feature = None
    
    for line in attrsfile:
        attr = line.rstrip("\r\n").split("=", 1)
        
        # everything before the "="
        cyto_info = attr[0].strip()
        feature = attr[1].strip().upper()
        fields = (cyto_info, "= True\n")
        
        if previous_feature and previous_feature == feature:
            fileout.write(" ".join(map(str, fields)))
        else:
            fileout = open("%s.%s.eda" % (uniqid, feature), 'w')
            fileout.write("feature%s\n" % feature)
            fileout.write(" ".join(map(str, fields)))
            previous_feature = feature


if __name__ == "__main__":
    main()