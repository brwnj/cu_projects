#! /usr/bin/env bash
SAMPLES=(2_ACATCG 2_CACTGT 2_CGTGAT 2_GCCTAA 2_TGGTCA           # 5
            3_ACATCG 3_CACTGT 3_CGTGAT 3_GCCTAA 3_TGGTCA)       # 10
PI=bennett
RESULTS=$HOME/projects/bennett/results/common
READS=$HOME/projects/bennett/data/common
BIN=$HOME/devel/repertoire/repertoire
R1PRIMERS=$READS/r1_primers.fasta
R2PRIMERS=$READS/r2_primers.fasta
