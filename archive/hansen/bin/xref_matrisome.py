#!/usr/bin/env python
# encoding: utf-8
"""
From fully annotated mouse matrisome, add rat protein reference, percent
ortholog, and annotations for the rat match.
"""
import argparse
import sys
from toolshed import reader, nopen
from collections import defaultdict
from itertools import izip

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def parse_uniprot_flat(uniprot):
    uid = ""
    records = defaultdict(lambda : defaultdict(list))
    for i, line in enumerate(nopen(uniprot)):
        line = line.rstrip("\r\n").split("   ", 1)
        code = line[0]
        
        try:
            attrs = line[1]
        except IndexError:
            # record delimiter
            uid = ""
            continue
            
        # if i > 100000:
        #     return records
        #     sys.exit()
        
        if code == 'ID':
            uid = attrs.split()[0]
        
        if uid == "":
            sys.stderr.write(">> Parsing error. UID not set.\n")
            sys.exit(1)
            
        if code == 'AC':
            uniprotid = attrs.rstrip(';').split("; ")
            [records[uid]['uniprotid'].append(u) for u in uniprotid]
        
        if code == 'DR':
            attrs = attrs.rstrip('.').split('; ')
            attrtype = attrs[0].lower()
            if attrtype == 'interpro':
                records[uid][attrtype].append(attrs[1])
            if attrtype == 'geneid':
                records[uid][attrtype].append(attrs[1])
            if attrtype == 'refseq':
                # RSN
                records[uid]['refseqn'].append(attrs[2])
                # RSP
                records[uid]['refseqp'].append(attrs[1])
            if attrtype == 'ensembl':
                records[uid]['ensemblt'].append(attrs[1])
                records[uid]['ensemblp'].append(attrs[2])
                records[uid]['ensemblg'].append(attrs[3])
        
        if code == 'GN':
            if attrs.find('Name=') > -1:
                attrs = attrs.rstrip(';').split('; ')
                records[uid]['genename'].append(attrs[0].split('=')[1])
        
        if code == 'DE':
            if attrs.find('Full=') > -1:
                attrs = attrs.rstrip(';').split('Full=')
                records[uid]['description'].append(attrs[1])
    return records


def get_args():
    """Process arguments."""
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('matrisome', help='matrisome file with header')
    p.add_argument('xref', help='cross-reference file from inparanoid')
    p.add_argument('xref_col', help='column name in matrisome xlxs dump')
    p.add_argument('uniprot', help='uniprot_sprot flat file')
    p.add_argument('-v', '--verbose', action='store_true', 
                    help='maximum verbosity')
    args = p.parse_args()
    return args


def get_xref(xref):
    """Makes `xd` from supplied table of mouse and rat orthologs in many-to-
    many relationship
    
    returns defaultdict(defaultdict(list))
    """
    xd = defaultdict(lambda : defaultdict(list))
    for e in reader(xref, header="uid bits species score name percent".split()):
        if e['species'].startswith("R.norv"):
            xd[e['uid']]['ratnames'].append(e['name'])
            xd[e['uid']]['ratscores'].append(e['score'])
        else:
            xd[e['uid']]['orthonames'].append(e['name'])
            xd[e['uid']]['orthoscores'].append(e['score'])
    return xd


def main():
    args = get_args()
    
    if args.verbose:
        sys.stderr.write(">> building gene orthology cross-reference...\n")
    xref = get_xref(args.xref)
    
    if args.verbose:
        sys.stderr.write(">> building uniprot library...\n")
    uniprot = parse_uniprot_flat(args.uniprot)
    
    if args.verbose:
        sys.stderr.write(">> annotating matrisome...\n")
    
    header = nopen(args.matrisome).readline().rstrip("\r\n").split("\t")
    headerext = ['r_ENSRNOP', 'r_score', 'r_geneid', 'r_gene_description', \
                    'r_uniprot', 'r_interpro', 'r_refseqn', 'r_refseqp', \
                    'r_ensg', 'r_enst', 'r_ensp']
    header.extend(headerext)
    print "\t".join(h for h in header)
    
    for entry in reader(args.matrisome):
        
        # reset vars
        for h in headerext:
            entry[h] = ""
        
        # handle multiple entries delimited by ":"
        for entryname in entry[args.xref_col].split(":"):
            
            # looping over entire defaultdict each time
            for uid, ddict in xref.iteritems():
                
                # find a matching ortholog
                for orthoname in ddict['orthonames']:
                    if orthoname == entryname:
                        
                        # use the uid to get the rat names and scores
                        for ratname, ratscore in izip(xref[uid]['ratnames'], xref[uid]['ratscores']):
                            # print ratname
                            entry['r_ENSRNOP'] += "%s:" % ratname
                            entry['r_score'] += "%s:" % ratscore
                            
                            # for each rat ENSP, add the corresponding annotation(s)
                            for uniqueid, uniprot_entry in uniprot.iteritems():
                                for ensemblname in uniprot_entry['ensemblp']:
                                    if ensemblname == ratname:
                                        #print all of the info for this uid
                                        entry['r_geneid'] += ':'.join(t for t in uniprot[uniqueid]['geneid']) + ":"
                                        entry['r_gene_description'] += ':'.join(t for t in uniprot[uniqueid]['description']) + ":"
                                        entry['r_uniprot'] += ':'.join(t for t in uniprot[uniqueid]['uniprotid']) + ":"
                                        entry['r_interpro'] += ':'.join(t for t in uniprot[uniqueid]['interpro']) + ":"
                                        entry['r_refseqn'] += ':'.join(t for t in uniprot[uniqueid]['refseqn']) + ":"
                                        entry['r_refseqp'] += ':'.join(t for t in uniprot[uniqueid]['refseqp']) + ":"
                                        entry['r_ensg'] += ':'.join(t for t in uniprot[uniqueid]['ensemblg']) + ":"
                                        entry['r_enst'] += ':'.join(t for t in uniprot[uniqueid]['ensemblt']) + ":"
                                        entry['r_ensp'] += ':'.join(t for t in uniprot[uniqueid]['ensemblp']) + ":"
                                        
        print "\t".join(entry[h].rstrip(":") for h in header)


if __name__ == "__main__":
    main()