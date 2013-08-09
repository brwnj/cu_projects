#!/usr/bin/env bash
#BSUB -J sra
#BSUB -e sra.%J.err
#BSUB -o sra.%J.out
#BSUB -q normal
#BSUB -R "select[mem>1] rusage[mem=1] span[hosts=1]"
#BSUB -n 1
#BSUB -P davidson

set -o nounset -o pipefail -o errexit -x

wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP010%2FSRP010483/SRR401815/SRR401815.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP010%2FSRP010483/SRR401816/SRR401816.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP010%2FSRP010483/SRR401817/SRR401817.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP010%2FSRP010483/SRR401818/SRR401818.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP010%2FSRP010483/SRR401819/SRR401819.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP010%2FSRP010483/SRR401820/SRR401820.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP010%2FSRP010483/SRR401821/SRR401821.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP010%2FSRP010483/SRR401822/SRR401822.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP010%2FSRP010483/SRR401823/SRR401823.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP010%2FSRP010483/SRR401824/SRR401824.sra
