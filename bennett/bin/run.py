#!/usr/bin/env python
# encoding: utf-8
"""
transfer miseq files to project/data/common
"""
import os
import urllib2
import pandas as pd
from BeautifulSoup import BeautifulSoup as bs

data = "/vol1/home/brownj/projects/bennett/data/common"
meta_tsv = "/vol1/home/brownj/projects/bennett/data/common/metadata.tsv"
meta = pd.read_table(meta_tsv, index_col=['patient', 'barcode'])

def find_files(barcode, soup):
    links = []
    for link in soup.find_all("a"):
        l = link.get("href")
        if barcode in l and l.endswith(".gz"):
            links.append[l]
    return links

def download(filename, url):
    f = urllib2.urlopen(url)
    with open(filename, 'wb') as fq:
        fq.write(f.read())

htmlsoup = {}
patients = set(meta.index.get_level_values('patient'))
for patient in patients:
    url = meta.ix[patient]['url'][0]
    response = urllib2.urlopen(url)
    html = response.read()
    htmlsoup[patient] = bs(html)

for index, row in meta.iterrows():
    patient, barcode = index
    file_links = find_files(barcode, htmlsoup[patient])
    index_read = "{data}/{patient}_{cell_type}_I1.fastq.gz".format(data=data, patient=patient, cell_type=row['cell_type'])
    r1 = "{data}/{patient}_{cell_type}_R1.fastq.gz".format(data=data, patient=patient, cell_type=row['cell_type'])
    r2 = "{data}/{patient}_{cell_type}_R2.fastq.gz".format(data=data, patient=patient, cell_type=row['cell_type'])
    # index read or I1
    if not os.path.exists(index_read):
        file_link = [f for f in file_links if "I1" in f][0]
        download(index_read, file_link)
    # R1
    if not os.path.exists(r1):
        file_link = [f for f in file_links if "R1" in f][0]
        download(r1, file_link)
    # R2
    if not os.path.exists(r2):
        file_link = [f for f in file_links if "R2" in f][0]
        download(r2, file_link)
