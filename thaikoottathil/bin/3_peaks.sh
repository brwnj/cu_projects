#! /usr/bin/env bash
#BSUB -J macs.peaks
#BSUB -e peaks.%J.err
#BSUB -o peaks.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1

<<DOC
use macs to call peaks.

In this technique, chromatin from each sample which is either untreated/Control (CT) or IL-13 treated samples

I wanted another comparison between untreated/ control vs IL-13 treated cells. I want to know whether IL-13 treatment alters STAT6 or c-Jun binding to the SPLUNC1 (Chr 20) or Muc5AC promoter region (Chr 6)?

No antibody control (No ab control) can serve as background for all the samples.
Each samples has its own CT IgG control also.

sr #	ChIP	Treatment/donor #
1	Input	CT -NJ #21- 2h
3	STAT6	CT -NJ #21- 2h

9	Input	CT -NJ #34- 2h
10	CT IgG	CT -NJ #34- 2h
11	STAT6	CT -NJ #34- 2h
12	c-Jun	CT -NJ #34- 2h

21	Input	CT 24h -NJ #34
22	STAT6	CT 24h -NJ #34
23	c-Jun	CT 24h -NJ #34

5	Input	IL-13 2h -NJ #21
7	STAT6	IL-13 2h -NJ #21

13	Input	IL-13 2h -NJ #34
14	CT IgG	IL-13 2h -NJ #34
15	STAT6	IL-13 2h -NJ #34
16	c-Jun	IL-13 2h -NJ #34

17	Input	IL-13 24h -NJ #34
18	STAT6	IL-13 24h -NJ #34
24	c-Jun	IL-13 24h -NJ #34

19	RNA pol II	RNA pol II ab (positive CT)
20	no antibody	No ab control  (Negative CT)

9_CTGATC_L006_R1_001
10_AAGCTA_L006_R1_001
11_GTAGCC_L006_R1_001
12_TACAAG_L006_R1_001
13_TTGACT_L006_R1_001
14_GGAACT_L006_R1_001
15_TGACAT_L006_R1_001
16_GGACGG_L006_R1_001
17_GCGGAC_L006_R1_001
18_TTTCAC_L006_R1_001
24_ATTGGC_L006_R1_001
19_GGCCAC_L006_R1_001
1_CGTGAT_L006_R1_001
20_CGAAAC_L006_R1_001
21_TGGTCA_L006_R1_001
22_TCAAGT_L006_R1_001
23_ACATCG_L006_R1_001
3_GCCTAA_L006_R1_001
5_CACTGT_L006_R1_001
7_GATCTG_L006_R1_001
DOC

set -o nounset -o pipefail -o errexit -x

PROJECT=/vol1/home/brownj/projects/thaikoottathil/results/common
# CONTROL=$PROJECT/20_CGAAAC_L006_R1_001/20_CGAAAC_L006_R1_001.bam

# for f in $PROJECT/*; do
#     if [ $f != 20_CGAAAC_L006_R1_001 ]; then
#         NAME=$(basename $f)
#         macs14 --treatment $PROJECT/$NAME/$NAME.bam \
#             --control $CONTROL \
#             --name $NAME \
#             --format BAM \
#             --gsize hs \
#             --wig \
#             --single-profile \
#             --call-subpeaks
#     fi
# done

macs14 --treatment $PROJECT/3_GCCTAA_L006_R1_001/3_GCCTAA_L006_R1_001.bam \
    --control $PROJECT/1_CGTGAT_L006_R1_001/1_CGTGAT_L006_R1_001.bam \
    --name 3_GCCTAA_L006_R1_001 \
    --format BAM \
    --gsize hs \
    --wig \
    --single-profile \
    --call-subpeaks

macs14 --treatment $PROJECT/7_GATCTG_L006_R1_001/7_GATCTG_L006_R1_001.bam \
    --control $PROJECT/5_CACTGT_L006_R1_001/5_CACTGT_L006_R1_001.bam \
    --name 7_GATCTG_L006_R1_001 \
    --format BAM \
    --gsize hs \
    --wig \
    --single-profile \
    --call-subpeaks

macs14 --treatment $PROJECT/10_AAGCTA_L006_R1_001/10_AAGCTA_L006_R1_001.bam \
    --control $PROJECT/9_CTGATC_L006_R1_001/9_CTGATC_L006_R1_001.bam \
    --name 10_AAGCTA_L006_R1_001 \
    --format BAM \
    --gsize hs \
    --wig \
    --single-profile \
    --call-subpeaks

macs14 --treatment $PROJECT/11_GTAGCC_L006_R1_001/11_GTAGCC_L006_R1_001.bam \
    --control $PROJECT/9_CTGATC_L006_R1_001/9_CTGATC_L006_R1_001.bam \
    --name 11_GTAGCC_L006_R1_001 \
    --format BAM \
    --gsize hs \
    --wig \
    --single-profile \
    --call-subpeaks

macs14 --treatment $PROJECT/12_TACAAG_L006_R1_001/12_TACAAG_L006_R1_001.bam \
    --control $PROJECT/9_CTGATC_L006_R1_001/9_CTGATC_L006_R1_001.bam \
    --name 12_TACAAG_L006_R1_001 \
    --format BAM \
    --gsize hs \
    --wig \
    --single-profile \
    --call-subpeaks

macs14 --treatment $PROJECT/14_GGAACT_L006_R1_001/14_GGAACT_L006_R1_001.bam \
    --control $PROJECT/13_TTGACT_L006_R1_001/13_TTGACT_L006_R1_001.bam \
    --name 14_GGAACT_L006_R1_001 \
    --format BAM \
    --gsize hs \
    --wig \
    --single-profile \
    --call-subpeaks

macs14 --treatment $PROJECT/15_TGACAT_L006_R1_001/15_TGACAT_L006_R1_001.bam \
    --control $PROJECT/13_TTGACT_L006_R1_001/13_TTGACT_L006_R1_001.bam \
    --name 15_TGACAT_L006_R1_001 \
    --format BAM \
    --gsize hs \
    --wig \
    --single-profile \
    --call-subpeaks

macs14 --treatment $PROJECT/16_GGACGG_L006_R1_001/16_GGACGG_L006_R1_001.bam \
    --control $PROJECT/13_TTGACT_L006_R1_001/13_TTGACT_L006_R1_001.bam \
    --name 16_GGACGG_L006_R1_001 \
    --format BAM \
    --gsize hs \
    --wig \
    --single-profile \
    --call-subpeaks

macs14 --treatment $PROJECT/18_TTTCAC_L006_R1_001/18_TTTCAC_L006_R1_001.bam \
    --control $PROJECT/17_GCGGAC_L006_R1_001/17_GCGGAC_L006_R1_001.bam \
    --name 18_TTTCAC_L006_R1_001 \
    --format BAM \
    --gsize hs \
    --wig \
    --single-profile \
    --call-subpeaks

macs14 --treatment $PROJECT/24_ATTGGC_L006_R1_001/24_ATTGGC_L006_R1_001.bam \
    --control $PROJECT/17_GCGGAC_L006_R1_001/17_GCGGAC_L006_R1_001.bam \
    --name 24_ATTGGC_L006_R1_001 \
    --format BAM \
    --gsize hs \
    --wig \
    --single-profile \
    --call-subpeaks

macs14 --treatment $PROJECT/22_TCAAGT_L006_R1_001/22_TCAAGT_L006_R1_001.bam \
    --control $PROJECT/21_TGGTCA_L006_R1_001/21_TGGTCA_L006_R1_001.bam \
    --name 22_TCAAGT_L006_R1_001 \
    --format BAM \
    --gsize hs \
    --wig \
    --single-profile \
    --call-subpeaks

macs14 --treatment $PROJECT/23_ACATCG_L006_R1_001/23_ACATCG_L006_R1_001.bam \
    --control $PROJECT/21_TGGTCA_L006_R1_001/21_TGGTCA_L006_R1_001.bam \
    --name 23_ACATCG_L006_R1_001 \
    --format BAM \
    --gsize hs \
    --wig \
    --single-profile \
    --call-subpeaks