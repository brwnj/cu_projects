#! /usr/bin/env bash
SAMPLES=(ON10-03_ACATCG ON10-03_ATTGGC ON10-03_CACTGT ON10-03_CGTGAT
         ON10-03_GATCTG ON10-03_GCCTAA ON10-03_TCAAGT ON10-03_TGGTCA
         MS13-02_ACATCG MS13-02_ATTGGC MS13-02_CACTGT MS13-02_CGTGAT
         MS13-02_GATCTG MS13-02_GCCTAA MS13-02_TCAAGT MS13-02_TGGTCA
         ON07-05_ACATCG ON07-05_ATTGGC ON07-05_CACTGT ON07-05_CGTGAT
         ON07-05_GATCTG ON07-05_GCCTAA ON07-05_TCAAGT ON07-05_TGGTCA
         ON09-03_ACATCG ON09-03_CACTGT ON09-03_CGTGAT ON09-03_GCCTAA
         ON09-03_TGGTCA ON10-01_ACATCG ON10-01_ATTGGC ON10-01_CACTGT
         ON10-01_CGTGAT ON10-01_GCCTAA ON10-01_TGGTCA)

PI=bennett
RESULTS=$HOME/projects/bennett/results/common
READS=$HOME/projects/bennett/data/common
BIN=$HOME/devel/repertoire/repertoire
R1PRIMERS=$READS/r1_primers.fasta
R2PRIMERS=$READS/r2_primers.fasta
UMILENGTH=8
UMISEQ=NNNNNNNN
