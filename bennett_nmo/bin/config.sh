#!/usr/bin/env bash
set -o nounset

SAMPLES=(
control_ACATCG
control_ATTGGC
control_CACTGT
control_CGTGAT
control_GATCTG
control_GCCTAA
control_TCAAGT
control_TGGTCA
ON07-05A_gc1
ON07-05A_memorya
ON07-05A_memoryb
ON07-05A_naive
ON07-05A_plasmablasta
ON07-05A_plasmablastb
ON07-05B_memory
ON07-05B_plasmablasts
ON08-08_gc
ON08-08_memory
ON08-08_naive
ON08-08_plasmablast
ON09-03A_gc1
ON09-03A_gc2
ON09-03A_memory
ON09-03A_naive
ON09-03A_plasmablast
ON09-03B_gc1
ON09-03B_gc2
ON09-03B_memory
ON09-03B_naive
ON09-03B_plasmablast
ON10-01A_gc1
ON10-01A_gc2
ON10-01A_memory
ON10-01A_naive
ON10-01A_plasmablast
ON10-01B_gc1
ON10-01B_gc2
ON10-01B_memory
ON10-01B_naive
ON10-01B_plasmablast
ON10-03A_gc1
ON10-03A_gc2
ON10-03A_memory
ON10-03A_naive
ON10-03A_plasmablast
ON10-03B_gc1
ON10-03B_gc2
ON10-03B_memory
ON10-03B_naive
ON10-03B_plasmablast
ON11-04_gc
ON11-04_memory
ON11-04_naive
ON11-04_plasmablasts
TUM09-527_gc
TUM09-527_plasmablast
)

PI=bennett
RESULTS=$HOME/projects/bennett_nmo/results/common
DATA=$HOME/projects/bennett_nmo/data/common
R1PRIMERS=$DATA/r1_primers.fasta
if [[ -f $R1PRIMERS.gz ]]; then
	gunzip -f $R1PRIMERS.gz
fi
R2PRIMERS=$DATA/r2_primers.fasta
if [[ -f $R2PRIMERS.gz ]]; then
	gunzip -f $R2PRIMERS.gz
fi
R2PRIMEROFFSET=$DATA/r2_primers_offsets-reverse.tab
UMILENGTH=8

PRCONS=0.6
MAXDIV=0.1
MINQUAL=20
R1_MAXERR=0.2
R2_MAXERR=0.2
AP_MAXERR=0.2
ALPHA=0.1
FS_MISS=10
CS_MISS=10
NPROC=12
