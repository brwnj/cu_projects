#! /usr/bin/env bash
#BSUB -J MIRA.assembly
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q bigmem
#BSUB -R "span[hosts=1]"

<<DOC
hybrid assembly using MIRA 3.4.0

# handling 454 data
$ python sff_extract.py -l orig_data/linkers.fasta -i "insert_size:5000,insert_stdev:900" -o saillant orig_data/HA5F6SE01.sff 
$ python sff_extract.py -l orig_data/linkers.fasta -i "insert_size:8000,insert_stdev:1000" -a -o saillant orig_data/HA5F6SE02.sff
$ mv saillant.fastq lcampe_in.454.fastq
$ mv saillant.xml lcampe_traceinfo_in.454.xml

# handling solexa data
$ awk '{if(/^@HWI/){print $0"/1"}else{print $0}}' orig_data/6_10_11_RS_CGATGT_L003_R1_001.fastq > lcampe_in.solexa.fastq
$ awk '{if(/^@HWI/){print $0"/1"}else{print $0}}' orig_data/6_10_11_RS_CGATGT_L003_R2_001.fastq >> lcampe_in.solexa.fastq

# read names are recommended to be shorter than 40 characters
$ echo "awk '{if(/^@HWI/){split(\$1,name,\":\");split(\$2,rnum,\"/\");print name[1]\":\"name[4]\":\"name[5]\":\"name[6]\":\"name[7]\"/\"rnum[2]}else{print \$0}}' lcampe_in.solexa.fastq > lcampe_in.solexa.fastq.fixedreadnames" | bsub-ez fix.read.names
DOC

set -o nounset -o pipefail -o errexit -x

mira --project=lcampe --job=denovo,genome,accurate,454,solexa -MI:sonfs=no SOLEXA_SETTINGS -GE:tismin=150:tismax=750 >&log_assembly.txt