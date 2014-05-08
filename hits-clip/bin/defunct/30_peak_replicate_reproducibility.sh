#! /usr/bin/env bash
#BSUB -J counting_peaks
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q normal

<<DOC
count the number of peaks while increasing stringency across cases across replicates
DOC

echo 'MDA231'
echo -n '1/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK24/PK24.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK42/PK42.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK54/PK54.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK24/PK24.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK42/PK42.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK54/PK54.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '2/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK24/PK24.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK42/PK42.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK54/PK54.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK24/PK24.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK42/PK42.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK54/PK54.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '3/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 3 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK24/PK24.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK42/PK42.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK54/PK54.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 3 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK24/PK24.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK42/PK42.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK54/PK54.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo ''
echo 'BT474'
echo -n '1/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK21/PK21.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK41/PK41.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK52/PK52.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK21/PK21.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK41/PK41.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK52/PK52.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '2/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK21/PK21.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK41/PK41.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK52/PK52.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK21/PK21.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK41/PK41.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK52/PK52.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '3/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 3 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK21/PK21.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK41/PK41.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK52/PK52.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 3 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK21/PK21.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK41/PK41.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK52/PK52.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo ''
echo 'BT474estr'
echo -n '1/2: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK22/PK22.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK53/PK53.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK22/PK22.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK53/PK53.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '2/2: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK22/PK22.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK53/PK53.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK22/PK22.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK53/PK53.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo ''
echo 'HS27A'
echo -n '1/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP1/MP1.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP21/MP21.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP35/MP35.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP1/MP1.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP21/MP21.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP35/MP35.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '2/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP1/MP1.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP21/MP21.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP35/MP35.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP1/MP1.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP21/MP21.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP35/MP35.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '3/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 3 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP1/MP1.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP21/MP21.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP35/MP35.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 3 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP1/MP1.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP21/MP21.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP35/MP35.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo ''
echo 'HS5'
echo -n '1/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP2/MP2.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP20/MP20.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP34/MP34.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP2/MP2.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP20/MP20.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP34/MP34.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '2/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP2/MP2.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP20/MP20.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP34/MP34.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP2/MP2.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP20/MP20.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP34/MP34.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '3/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 3 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP2/MP2.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP20/MP20.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP34/MP34.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 3 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP2/MP2.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP20/MP20.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP34/MP34.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo ''
echo 'hMSC'
echo -n '1/5: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP36/MP36.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.ACTG/MP43.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.TCGA/MP43.TCGA.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.ACTG/MP44.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.TCGA/MP44.TCGA.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP36/MP36.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.ACTG/MP43.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.TCGA/MP43.TCGA.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.ACTG/MP44.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.TCGA/MP44.TCGA.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '2/5: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP36/MP36.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.ACTG/MP43.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.TCGA/MP43.TCGA.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.ACTG/MP44.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.TCGA/MP44.TCGA.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP36/MP36.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.ACTG/MP43.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.TCGA/MP43.TCGA.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.ACTG/MP44.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.TCGA/MP44.TCGA.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '3/5: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 3 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP36/MP36.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.ACTG/MP43.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.TCGA/MP43.TCGA.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.ACTG/MP44.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.TCGA/MP44.TCGA.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 3 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP36/MP36.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.ACTG/MP43.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.TCGA/MP43.TCGA.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.ACTG/MP44.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.TCGA/MP44.TCGA.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '4/5: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 4 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP36/MP36.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.ACTG/MP43.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.TCGA/MP43.TCGA.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.ACTG/MP44.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.TCGA/MP44.TCGA.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 4 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP36/MP36.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.ACTG/MP43.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.TCGA/MP43.TCGA.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.ACTG/MP44.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.TCGA/MP44.TCGA.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '5/5: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 5 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP36/MP36.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.ACTG/MP43.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.TCGA/MP43.TCGA.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.ACTG/MP44.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.TCGA/MP44.TCGA.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 5 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP36/MP36.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.ACTG/MP43.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP43.TCGA/MP43.TCGA.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.ACTG/MP44.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP44.TCGA/MP44.TCGA.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo ''
echo 'HUVEC'
echo -n '1/2: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP24/MP24.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP38/MP38.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP24/MP24.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP38/MP38.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '2/2: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP24/MP24.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP38/MP38.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP24/MP24.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP38/MP38.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo ''
echo 'MCF7estr'
echo -n '1/2: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK12/PK12.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK32/PK32.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK12/PK12.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK32/PK32.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '2/2: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK12/PK12.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK32/PK32.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK12/PK12.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK32/PK32.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo ''
echo 'BT474herc'
echo 'no replicates'
echo ''
echo 'BMEC'
echo -n '1/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP42.ACTG/MP42.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP45.ACTG/MP45.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP45.TCGA/MP45.TCGA.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP42.ACTG/MP42.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP45.ACTG/MP45.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP45.TCGA/MP45.TCGA.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '2/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP42.ACTG/MP42.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP45.ACTG/MP45.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP45.TCGA/MP45.TCGA.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP42.ACTG/MP42.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP45.ACTG/MP45.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP45.TCGA/MP45.TCGA.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '3/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 3 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP42.ACTG/MP42.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP45.ACTG/MP45.ACTG.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP45.TCGA/MP45.TCGA.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 3 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP42.ACTG/MP42.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP45.ACTG/MP45.ACTG.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/MP45.TCGA/MP45.TCGA.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo ''
echo 'MCF7'
echo -n '1/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK11/PK11.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK31/PK31.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK51/PK51.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 1 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK11/PK11.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK31/PK31.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK51/PK51.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '2/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK11/PK11.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK31/PK31.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK51/PK51.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 2 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK11/PK11.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK31/PK31.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK51/PK51.pos.peaks.rmd.qv.001.bed.gz | wc -l))
echo -n '3/3: '
echo $(expr $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 3 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK11/PK11.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK31/PK31.neg.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK51/PK51.neg.peaks.rmd.qv.001.bed.gz | wc -l) \
    + \
    $(python ~/devel/peaktools/peaktools/combine_replicates.py -m 3 \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK11/PK11.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK31/PK31.pos.peaks.rmd.qv.001.bed.gz \
    /vol1/home/brownj/projects/hits-clip/results/common/samples/PK51/PK51.pos.peaks.rmd.qv.001.bed.gz | wc -l))