#!/usr/bin/env python
# encoding: utf-8
"""
parses chippeakanno results, adding chr?:start-stop column and removing "names"
column.
"""
import sys
from toolshed import reader

def main(args):
    header = "chrom start stop browserhelp width strand feature start_position \
                stop_position insideFeature distancetoFeature shortestDistance \
                fromOverlappingOrNearest symbol".split()
    print "\t".join(header)
    for t in reader(args.result):
        fields = ["chr%s" % t["space"], t["start"], t["end"],
                    "chr%s:%s-%s" % (t['space'], t['start'], t['end']),
                    t['width'], t['strand'], t['feature'],
                    t['start_position'], t['end_position'], t['insideFeature'], 
                    t['distancetoFeature'], t['shortestDistance'], 
                    t['fromOverlappingOrNearest'], t['symbol']]
        print "\t".join(map(str, fields))

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('result', 
                help='final result of chippeakanno after adding gene symbol')
    args = p.parse_args()
    main(args)