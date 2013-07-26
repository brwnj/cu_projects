#! /usr/bin/env bash
SAMPLES=(2_ACATCG 2_CACTGT 2_CGTGAT 2_GCCTAA 2_TGGTCA           # 5
            3_ACATCG 3_CACTGT 3_CGTGAT 3_GCCTAA 3_TGGTCA        # 10
            4_ACATCG 4_ATTGGC 4_CACTGT 4_CGTGAT 4_GATCTG        # 15
            4_GCCTAA 4_TCAAGT 4_TGGTCA)                         # 18
PI=bennett
RESULTS=$HOME/projects/bennett/results/common
READS=$HOME/projects/bennett/data/common
BIN=$HOME/devel/repertoire/repertoire
R1PRIMERS=$READS/r1_primers.fasta
R2PRIMERS=$READS/r2_primers.fasta
UMILENGTH=8
UMISEQ=NNNNNNNN
