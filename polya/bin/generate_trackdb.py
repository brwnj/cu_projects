#!/usr/bin/env python
# encoding: utf-8
"""
Create trackdb from list of folder contents.
"""
import os
import re
import sys
import argparse
from string import Template
from toolshed import reader
from collections import defaultdict

COMPOSITE_TRACK_DEF = ("track PolyA\n"
                        "compositeTrack on\n"
                        "shortLabel Sites and Coverage\n"
                        "longLabel Sites and Coverage tracks\n"
                        "subGroup1 view Views PKS=Classified_Sites COV=Coverage SHT=Shifts STS=Sites\n"
                        "subGroup2 stype SampleType UNK=Unknown CD71PCD235P=CD71+/CD235a+ CD71PCD235N=CD71+/CD235a- CD71NCD235P=CD71-/CD235a+ CD34=CD34 TS=Tumor NBT=NormalTissue\n"
                        "subGroup3 strand Strand POS=Positive NEG=Negative\n"
                        "subGroup4 inv Investigator MP=Pillai PK=Kabos\n"
                        "subGroup5 day Day 1=1 4=4 7=7 14=14\n"
                        "sortOrder view=-\n"
                        "type bed 6 .\n")

CLASSIFIED_PEAKS_VIEW = ("  track viewPeaks\n"
                            "   parent PolyA\n"
                            "   shortLabel Classified Sites\n"
                            "   view PKS\n"
                            "   visibility pack\n"
                            "   type bigBed 6 .\n")

COVERAGE_VIEW = ("  track viewCoverage\n"
                    "   parent PolyA\n"
                    "   shortLabel Coverage\n"
                    "   view COV\n"
                    "   visibility full\n"
                    "   maxHeightPixels 50:20:15\n"
                    "   type bigWig\n"
                    "   autoScale on\n"
                    "   alwaysZero on\n")

SITES_VIEW = """
    track viewSites
    parent PolyA
    shortLabel Sites
    view STS
    visibility pack
    type bigBed 6 .

        track mp.c13.sites
        parent viewSites
        subGroups view=STS inv=MP
        shortLabel MPSites13
        longLabel MP sites classes 1 and 3
        bigDataUrl MP.sites.c13.bb
        color 215,48,39
        type bigBed 6 .

        track mp.c1234.sites
        parent viewSites
        subGroups view=STS inv=MP
        shortLabel MPSites1234
        longLabel MP sites classes 1, 2, 3, and 4
        bigDataUrl MP.sites.c1234.bb
        color 215,48,39
        type bigBed 6 .

        track mp.wholegene.sites
        parent viewSites
        subGroups view=STS inv=MP
        shortLabel MPSitesWG
        longLabel MP sites classes 1, 2, 3, and 4 - Whole Gene
        bigDataUrl MP.sites.wholegene.bb
        color 189,0,38
        type bigBed 6 .

        track pk.c13.sites
        parent viewSites
        subGroups view=STS inv=PK
        shortLabel PKSites13
        longLabel PK sites classes 1 and 3
        bigDataUrl PK.sites.c13.bb
        color 69,117,180
        type bigBed 6 .

        track pk.c1234.sites
        parent viewSites
        subGroups view=STS inv=PK
        shortLabel PKSites1234
        longLabel PK sites classes 1, 2, 3, and 4
        bigDataUrl PK.sites.c1234.bb
        color 69,117,180
        type bigBed 6 .

        track pk.wholegene.sites
        parent viewSites
        subGroups view=STS inv=PK
        shortLabel PKSitesWG
        longLabel PK sites classes 1, 2, 3, and 4 - Whole Gene
        bigDataUrl PK.sites.wholegene.bb
        color 0,109,44
        type bigBed 6 .
"""

# these will presumably have different sample subtypes
PK_SHIFTS = ("track pkshifts\n"
                "compositeTrack on\n"
                "configurable on\n"
                "shortLabel PK Shifts\n"
                "longLabel Kabos: Observed DEXSeq shifts\n"
                "subGroup1 stype SampleType UNK=Unknown CD71P=CD71pos CD71N=CD71neg CD235AP=CD235Apos CD235AN=CD235Aneg TS=Tumor NBT=NormalTissue\n"
                "subGroup2 strand Strand POS=Positive NEG=Negative\n"
                "type bed 12 .\n")

MP_SHIFTS = ("track mpshifts\n"
                "compositeTrack on\n"
                "configurable on\n"
                "shortLabel MP Shifts\n"
                "longLabel Pillai: Observed DEXSeq shifts\n"
                "subGroup1 stype SampleType UNK=Unknown CD71P=CD71pos CD71N=CD71neg CD235AP=CD235Apos CD235AN=CD235Aneg TS=Tumor NBT=NormalTissue\n"
                "subGroup2 strand Strand POS=Positive NEG=Negative\n"
                "type bed 12 .\n")

SITE_TEMPLATE = Template("        track $tname\n"
                            "        parent viewPeaks\n"
                            "        subGroups stype=$stype view=PKS inv=$inv\n"
                            "        shortLabel $tname\n"
                            "        longLabel $tname\n"
                            "        bigDataUrl $filename\n"
                            "        color $color\n"
                            "        type bigBed 6 .\n")

COVERAGE_TEMPLATE = Template("        track $tname\n"
                                "        parent viewCoverage\n"
                                "        subGroups stype=$stype view=COV inv=$inv strand=$strand\n"
                                "        shortLabel $tname\n"
                                "        longLabel $tname\n"
                                "        bigDataUrl $filename\n"
                                "        color $color\n"
                                "        type bigWig\n")

SHIFTS_TEMPLATE = Template("    track $tname\n"
                                "    parent $parent\n"
                                "    subGroups stype=$stype strand=$strand\n"
                                "    shortLabel $tname\n"
                                "    longLabel $tname\n"
                                "    bigDataUrl $filename\n"
                                "    color 37,52,148\n"
                                "    exonArrows on\n"
                                "    exonArrowsDense on\n"
                                '    nextExonText "Next Shift"\n'
                                '    prevExonText "Previous Shift"\n'
                                "    thickDrawItem on\n"
                                "    type bigBed 12 .\n")

def gsample(fname):
    return re.split('\.|_', fname)[0]

def gtname(fname):
    return os.path.splitext(fname)[0]

def ginv(sname):
    return "PK" if "PK" in sname else "MP"

def gstrand(fname):
    return "NEG" if "neg" in fname else "POS"

def gshiftsparent(fname):
    return "pkshifts" if fname.startswith("PK") else "mpshifts"

def flipstrand(fname):
    # since reads are being mapped to the opposite strand they belong to,
    # we are flipping the names, eg. 'pos' to 'neg'.
    tname = os.path.splitext(fname)[0]
    if "neg" in tname:
        tname = tname.replace("neg", "pos")
    elif "pos" in tname:
        tname = tname.replace("pos", "neg")
    return tname

def main(folder, meta):
    filelist = os.listdir(folder)
    # group files by track type
    files = defaultdict(list)
    for f in filelist:
        if f.endswith("bw"):
            files['coverage'].append(f)
            continue
        if "_to_" in f:
            files['dexseq'].append(f)
            continue
        if "classified" in f:
            files['sites'].append(f)
            continue
    # metadata
    md = {}
    for l in reader(meta):
        md[l['alias']] = l
    # composite track with classified peaks, coverage, and merged sites
    print COMPOSITE_TRACK_DEF

    # classified peaks for peter and manoj
    # $results/*_peaks.classified.bed.gz >> folder/${sample}_peaks.classified.bb
    print CLASSIFIED_PEAKS_VIEW
    for s in files['sites']:
        sample = gsample(s)
        print SITE_TEMPLATE.substitute(tname=gtname(s),
                                    stype=md[sample]['stype'],
                                    inv=ginv(s),
                                    filename=s,
                                    color=md[sample]['color'])
    print COVERAGE_VIEW
    # the bigwigs
    for s in files['coverage']:
        sample = gsample(s)
        tname = flipstrand(s)
        print COVERAGE_TEMPLATE.substitute(tname=tname,
                                    stype=md[sample]['stype'],
                                    inv=ginv(s),
                                    strand=gstrand(tname),
                                    filename=s,
                                    color=md[sample]['color'])
    print SITES_VIEW
    # should not change often, but could automate eventually
    print PK_SHIFTS
    # after classifying shifts
    for s in files['dexseq']:
        if not s.startswith("PK"): continue
        tname = flipstrand(s)
        print SHIFTS_TEMPLATE.substitute(tname=tname,
                                    parent=gshiftsparent(s),
                                    stype="UNK",
                                    strand=gstrand(tname),
                                    filename=s)
    print MP_SHIFTS
    for s in files['dexseq']:
        if not s.startswith("MP"): continue
        tname = flipstrand(s)
        print SHIFTS_TEMPLATE.substitute(tname=tname,
                                    parent=gshiftsparent(s),
                                    stype="UNK",
                                    strand=gstrand(tname),
                                    filename=s)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("folder", help="files to convert to track definitions")
    p.add_argument("meta", help="google docs with sample info as tsv")
    args = p.parse_args()
    main(args.folder, args.meta)
