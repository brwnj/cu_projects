#!/usr/bin/env bash
#BSUB -J reprocessing
#BSUB -e reprocessing.%J.err
#BSUB -o reprocessing.%J.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P hits-clip

<<DOC
The CLIP samples to filter are:

PK11, PK31 and PK15 (MCF-7)
PK21, PK41 and PK52 (BT474)
P24, PK42 and PK54 (MDA-231)

I think we should do this slightly differently than with the newer samples:

    I'd take each aligned read, extend it 25bp downstream (rather than 100bp),
    and look for the presence of any of the following sequences in the extended
    alignment:

        GTGTCA
        GTGCCA
        GTGTCT
        GTCTCA
        GTGACA

    Discard any reads with matches in the extended alignment, and then call
    peaks as before (including combined peaks for each cell line).
DOC

set -o nounset -o pipefail -o errexit -x

SAMPLES=(PK11 PK31 PK51 PK21 PK41 PK52 PK24 PK42 PK54)

fasta=$HOME/ref/hg19/hg19.fa
fasta_dir=$HOME/ref/hg19/fa_per_chr
sizes=$HOME/ref/hg19/hg19.sizes
results=/vol1/home/brownj/projects/hits-clip/results/20140317

groups=(MCF7 BT474 MDA231)
declare -A replicates
replicates=([MCF7]="PK11 PK31 PK51"
            [BT474]="PK21 PK41 PK52"
            [MDA231]="PK24 PK42 PK54")
declare -A colors
colors=([PK11]="0,109,44"
        [PK31]="35,139,69"
        [PK51]="65,174,118"
        [MCF7]="0,68,27"
        [PK21]="8,81,156"
        [PK41]="33,113,181"
        [PK52]="66,146,198"
        [BT474]="8,48,107"
        [PK24]="166,54,3"
        [PK42]="217,72,1"
        [PK54]="241,105,19"
        [MDA231]="127,39,4")
declare -A ucsc_group
ucsc_group=([PK11]=MCF7
            [PK31]=MCF7
            [PK51]=MCF7
            [PK21]=BT474
            [PK41]=BT474
            [PK52]=BT474
            [PK24]=MDA231
            [PK42]=MDA231
            [PK54]=MDA231
            [MDA231]=MDA231
            [MCF7]=MCF7
            [BT474]=BT474)

HUB=hub
GENOME=hg19
trackdb=$HUB/$GENOME/trackDb.txt

# bam1=/vol1/home/brownj/projects/hits-clip/results/common/$sample/$sample.bam
# bam2=/vol1/home/brownj/projects/hits-clip/results/common/$sample/$sample.rmd.bam
#
# outbam1=$results/${sample}_filtered.bam
# outbam2=$results/${sample}_filtered.rmd.bam
#
# python ~/projects/hits-clip/bin/scripts/filter_bam.py -m all -d GTGTCA -d GTGCCA -d GTGTCT -d GTCTCA -d GTGACA -b 25 $bam1 $outbam1 $fasta
# python ~/projects/hits-clip/bin/scripts/filter_bam.py -m all -d GTGTCA -d GTGCCA -d GTGTCT -d GTCTCA -d GTGACA -b 25 $bam2 $outbam2 $fasta

# build a new genomedata archive
# bam2gd.py -p hits-clip $sizes $fasta_dir *.bam

# call peaks
GENOMEDATA=/vol1/home/brownj/projects/hits-clip/results/20140317/genomedata

# peaks
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}
    non_unique_q=0.00001
    unique_q=0.01
    bam=${sample}_filtered.bam
    rmdbam=${sample}_filtered.rmd.bam

    # disaster.
    negbam=${sample}_filtered.neg.bam
    negNP=${negbam/.bam/_peaks.narrowPeak.gz}
    negbw=$HUB/$GENOME/${sample}_filtered_neg.bw
    negrmdbam=${sample}_filtered.rmd.neg.bam
    negrmdNP=${negrmdbam/.bam/_peaks.narrowPeak.gz}
    negrmdbw=$HUB/$GENOME/${sample}_filtered.rmd_neg.bw

    posbam=${sample}_filtered.pos.bam
    posNP=${posbam/.bam/_peaks.narrowPeak.gz}
    posbw=$HUB/$GENOME/${sample}_filtered_pos.bw
    posrmdbam=${sample}_filtered.rmd.pos.bam
    posrmdNP=${posrmdbam/.bam/_peaks.narrowPeak.gz}
    posrmdbw=$HUB/$GENOME/${sample}_filtered.rmd_pos.bw

    # stranded bams
    running_jobs=false
    if [[ ! -f $negbam ]]; then
        cmd="samtools view -hb -f 0x10 $bam > $negbam"
        bsub -J stranded_bam -o sb.%J.out -e sb.%J.err -P hits-clip -K $cmd &
        running_jobs=true
    fi
    if [[ ! -f $negrmdbam ]]; then
        cmd="samtools view -hb -f 0x10 $rmdbam > $negrmdbam"
        bsub -J stranded_bam -o sb.%J.out -e sb.%J.err -P hits-clip -K $cmd &
        running_jobs=true
    fi
    if [[ ! -f $posbam ]]; then
        cmd="samtools view -hb -F 0x10 $bam > $posbam"
        bsub -J stranded_bam -o sb.%J.out -e sb.%J.err -P hits-clip -K $cmd &
        running_jobs=true
    fi
    if [[ ! -f $posrmdbam ]]; then
        cmd="samtools view -hb -F 0x10 $rmdbam > $posrmdbam"
        bsub -J stranded_bam -o sb.%J.out -e sb.%J.err -P hits-clip -K $cmd &
        running_jobs=true
    fi
    if [[ running_jobs = true ]]; then
        wait
    fi

    # stranded bigwigs
    if [[ ! -f $negbw ]]; then
        t=$sample.tmp
        bedtools genomecov -strand - -bg -ibam $negbam | bedtools sort -i - > $t
        bedGraphToBigWig $t $sizes $negbw
        rm -f $t
    fi
    if [[ ! -f $negrmdbw ]]; then
        t=$sample.tmp
        bedtools genomecov -strand - -bg -ibam $negrmdbam | bedtools sort -i - > $t
        bedGraphToBigWig $t $sizes $negrmdbw
        rm -f $t
    fi
    if [[ ! -f $posbw ]]; then
        t=$sample.tmp
        bedtools genomecov -strand + -bg -ibam $posbam | bedtools sort -i - > $t
        bedGraphToBigWig $t $sizes $posbw
        rm -f $t
    fi
    if [[ ! -f $posrmdbw ]]; then
        t=$sample.tmp
        bedtools genomecov -strand + -bg -ibam $posrmdbam | bedtools sort -i - > $t
        bedGraphToBigWig $t $sizes $posrmdbw
        rm -f $t
    fi

    # call peaks on the stranded bams
    if [[ ! -f $negNP ]]; then
        cmd="macs2 callpeak -t $negbam -g hs -n ${sample}_filtered.neg --nomodel --extsize 20 -q $non_unique_q --keep-dup all"
        bsub -J peaks -o peaks.%J.out -e peaks.%J.err -P hits-clip -K $cmd &
    fi
    if [[ ! -f $negrmdNP ]]; then
        cmd="macs2 callpeak -t $negrmdbam -g hs -n ${sample}_filtered.rmd.neg --nomodel --extsize 20 -q $unique_q --keep-dup all"
        bsub -J peaks -o peaks.%J.out -e peaks.%J.err -P hits-clip -K $cmd &
    fi
    if [[ ! -f $posNP ]]; then
        cmd="macs2 callpeak -t $posbam -g hs -n ${sample}_filtered.pos --nomodel --extsize 20 -q $non_unique_q --keep-dup all"
        bsub -J peaks -o peaks.%J.out -e peaks.%J.err -P hits-clip -K $cmd &
    fi
    if [[ ! -f $posrmdNP ]]; then
        cmd="macs2 callpeak -t $posrmdbam -g hs -n ${sample}_filtered.rmd.pos --nomodel --extsize 20 -q $unique_q --keep-dup all"
        bsub -J peaks -o peaks.%J.out -e peaks.%J.err -P hits-clip -K $cmd &
    fi
done
echo ">> Waiting on Peak calls to finish"
wait

# clean up after MACS
if test -n "$(find . -maxdepth 1 -name '*.narrowPeak' -print -quit)"; then
    gzip -f *.narrowPeak
fi
if test -n "$(find . -maxdepth 1 -name '*.bed' -print -quit)"; then
    gzip -f *.bed
fi
if test -n "$(find . -maxdepth 1 -name '*.xls' -print -quit)"; then
    rm -f *.xls
fi

# merge peaks across replicates
for group in MCF7 BT474 MDA231; do
    for strand in pos neg; do
        tracks=""
        bedfiles=""
        combined=${group}_${strand}_rmd_peaks.bed.gz
        trimmed=$results/${group}_${strand}_rmd_trimmed_peaks.bed.gz

        for sample in ${replicates[$group]}; do
            # PK24_filtered.rmd.neg_peaks.narrowPeak.gz
            bedfiles="$bedfiles ${sample}_filtered.rmd.${strand}_peaks.narrowPeak.gz"
            tracks="$tracks -t ${sample}_filtered.rmd_${strand}"
        done

        # combined the replicates
        if [[ ! -f $combined ]]; then
            peaktools-combine-replicates --verbose $bedfiles | bedClip stdin $sizes ${combined/.gz}
            gzip -f ${combined/.gz}
        fi
        # submit trim job separately
        if [[ ! -f $trimmed ]]; then
            cmd="python ~/projects/hits-clip/bin/scripts/trim_peaks.py -v $tracks $combined $GENOMEDATA | bedtools sort -i - | gzip -c > $trimmed"
            bsub -J trim -o trim.%J.out -e trim.%J.err -P hits-clip -K $cmd &
        fi


        # also run bams containing duplicates
        tracks=""
        bedfiles=""
        combined=${group}_${strand}_peaks.bed.gz
        trimmed=${group}_${strand}_trimmed_peaks.bed.gz

        for sample in ${replicates[$group]}; do
            # PK24_filtered.neg_peaks.narrowPeak.gz
            bedfiles="$bedfiles ${sample}_filtered.${strand}_peaks.narrowPeak.gz"
            tracks="$tracks -t ${sample}_filtered_${strand}"
        done
        if [[ ! -f $combined ]]; then
            peaktools-combine-replicates --verbose $bedfiles | bedClip stdin $sizes ${combined/.gz}
            gzip -f ${combined/.gz}
        fi
        if [[ ! -f $trimmed ]]; then
            cmd="python ~/projects/hits-clip/bin/scripts/trim_peaks.py -v $tracks $combined $GENOMEDATA | bedtools sort -i - | gzip -c > $trimmed"
            bsub -J trim -o trim.%J.out -e trim.%J.err -P hits-clip -K $cmd &
        fi
    done
done
echo ">> Waiting on trim jobs to finish"
wait


if [[ ! -d $HUB/$GENOME ]]; then
    mkdir -p $HUB/$GENOME
fi

# genomes.txt
if [[ ! -f $HUB/genomes.txt ]]; then
    genomes=$HUB/genomes.txt
    echo "genome $GENOME" > $genomes
    echo "trackDb $GENOME/trackDb.txt" >> $genomes
fi

# hub.txt
if [[ ! -f $HUB/hub.txt ]]; then
    hub=$HUB/hub.txt
    echo "hub re_hitsclip" > $hub
    echo "shortLabel re-HITS-CLIP" >> $hub
    echo "longLabel re-HITS-CLIP" >> $hub
    echo "genomesFile genomes.txt" >> $hub
    echo "email brwnjm@gmail.com" >> $hub
fi

# output for coverage
cat <<coverage_track >$trackdb
track re_hitsclip_coverage
compositeTrack on
subGroup1 cline CellLine MDA231=MDA231 BT474=BT474 MCF7=MCF7
shortLabel Coverage
longLabel Coverage
maxHeightPixels 50:20:15
type bigWig
configurable on
autoScale on

coverage_track

for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}

    # already written into the hub directory
    posbw=${sample}_filtered_pos.bw
    negbw=${sample}_filtered_neg.bw
    posuniqbw=${sample}_filtered.rmd_pos.bw
    neguniqbw=${sample}_filtered.rmd_neg.bw

    color=${colors[$sample]}
    cat <<coverage_track >>$trackdb
    track ${posbw/.bw}
    bigDataUrl $posbw
    subGroups cline=${ucsc_group[$sample]}
    shortLabel $sample coverage POS
    longLabel $sample coverage positive (+) strand
    type bigWig
    parent re_hitsclip_coverage
    color $color

    track ${negbw/.bw}
    bigDataUrl $negbw
    subGroups cline=${ucsc_group[$sample]}
    shortLabel $sample coverage NEG
    longLabel $sample coverage negative (-) strand
    type bigWig
    parent re_hitsclip_coverage
    color $color

    track ${posuniqbw/.bw}
    bigDataUrl $posuniqbw
    subGroups cline=${ucsc_group[$sample]}
    shortLabel $sample unique coverage POS
    longLabel $sample unique coverage positive (+) strand
    type bigWig
    parent re_hitsclip_coverage
    color $color

    track ${neguniqbw/.bw}
    bigDataUrl $neguniqbw
    subGroups cline=${ucsc_group[$sample]}
    shortLabel $sample unique coverage NEG
    longLabel $sample unique coverage negative (-) strand
    type bigWig
    parent re_hitsclip_coverage
    color $color

coverage_track
done

# intervals
cat <<intervals_track >>$trackdb
track peak_intervals
compositeTrack on
shortLabel Peak Intervals
longLabel Peak Intervals
subGroup1 cline CellLine MDA231=MDA231 BT474=BT474 MCF7=MCF7
type bigBed 6

intervals_track

fields=bb_fields.as
cat <<bigbedfields >$fields
table hg19intervals
"NarrowPeak format"
(
string  chrom;		"Reference Chromosome"
uint    chromStart;	"Start Position"
uint    chromEnd;	"End Position"
string  name;		"Gene Name"
uint    score;		"Score"
char[1] strand;		"Strand: + or -"
float   signalValue;	"Measurement of Overall Enrichment"
float   pValue;	    "p-value (-log10)"
float  	qValue;	    "q-value (-log10)"
uint    peak;	    "Peak source (derived peak summit)"
)
bigbedfields

# individual samples
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}
    color=${colors[$sample]}

    # non-unique
    np_us=${sample}_filtered.narrowPeak.unsorted
    np=${sample}_filtered.narrowPeak
    bb=${np/.narrowPeak/.bb}
    awk -t -cbed '{if($5>1000){$5=1000}; $6="-"; print}' ${sample}_filtered.neg_peaks.narrowPeak.gz > $np_us
    awk -t -cbed '{if($5>1000){$5=1000}; $6="+"; print}' ${sample}_filtered.pos_peaks.narrowPeak.gz >> $np_us
    bedtools sort -i $np_us > $np
    bedToBigBed -type=bed6+4 -as=$fields $np $sizes $bb
    mv $bb $HUB/$GENOME
    rm -f $np_us
    gzip -f $np

    cat <<intervals_track >>$trackdb
        track ${bb/.bb}
        bigDataUrl $bb
        shortLabel ${bb/.bb} non-unique
        longLabel $sample: non-unique peaks
        subGroups cline=${ucsc_group[$sample]}
        type bigBed 6 +
        color $color
        parent peak_intervals

intervals_track

    # unique
    np_us=${sample}_filtered_rmd.narrowPeak.unsorted
    np=${sample}_filtered_rmd.narrowPeak
    bb=${np/.narrowPeak/.bb}
    awk -t -cbed '{if($5>1000){$5=1000}; $6="-"; print}' ${sample}_filtered.rmd.neg_peaks.narrowPeak.gz > $np_us
    awk -t -cbed '{if($5>1000){$5=1000}; $6="+"; print}' ${sample}_filtered.rmd.pos_peaks.narrowPeak.gz >> $np_us
    bedtools sort -i $np_us > $np
    bedToBigBed -type=bed6+4 -as=$fields $np $sizes $bb
    mv $bb $HUB/$GENOME
    rm -f $np_us
    gzip -f $np

    cat <<intervals_track >>$trackdb
        track ${bb/.bb}
        bigDataUrl $bb
        shortLabel ${bb/.bb} unique
        longLabel $sample: unique peaks
        subGroups cline=${ucsc_group[$sample]}
        type bigBed 6 +
        color $color
        parent peak_intervals

intervals_track

done

# combine stranded peaks into a single bigbed track
for group in BT474 MCF7 MDA231; do
    # unique
    awk -v group=$group -t -cbed '{print $1, $2, $3, group, "0", "+"}' ${group}_pos_rmd_trimmed_peaks.bed.gz >> ${group}_rmd_trimmed_peaks.unsorted
    awk -v group=$group -t -cbed '{print $1, $2, $3, group, "0", "-"}' ${group}_neg_rmd_trimmed_peaks.bed.gz >> ${group}_rmd_trimmed_peaks.unsorted
    bedtools sort -i ${group}_rmd_trimmed_peaks.unsorted | bedtools merge -s -i - | bedtools sort -i - | awk -v group=$group -t '{print $1, $2, $3, group, 0, $6}' > ${group}_rmd_trimmed_peaks.bed
    bedToBigBed -type=bed6 ${group}_rmd_trimmed_peaks.bed $sizes ${group}_rmd_trimmed_peaks.bb
    mv ${group}_rmd_trimmed_peaks.bb $HUB/$GENOME
    # non-unique
    awk -v group=$group -t -cbed '{print $1, $2, $3, group, "0", "+"}' ${group}_pos_trimmed_peaks.bed.gz >> ${group}_trimmed_peaks.unsorted
    awk -v group=$group -t -cbed '{print $1, $2, $3, group, "0", "-"}' ${group}_neg_trimmed_peaks.bed.gz >> ${group}_trimmed_peaks.unsorted
    bedtools sort -i ${group}_trimmed_peaks.unsorted | bedtools merge -s -i - | bedtools sort -i - | awk -v group=$group -t '{print $1, $2, $3, group, 0, $6}' > ${group}_trimmed_peaks.bed
    bedToBigBed -type=bed6 ${group}_trimmed_peaks.bed $sizes ${group}_trimmed_peaks.bb
    mv ${group}_trimmed_peaks.bb $HUB/$GENOME

    color=${colors[$group]}
    cat <<intervals_track >>$trackdb
        track ${group}_uintervals
        bigDataUrl ${group}_rmd_trimmed_peaks.bb
        shortLabel $group unique
        longLabel $group: Trimmed peak regions from unique reads
        subGroups cline=${ucsc_group[$group]}
        type bigBed 6
        color $color
        parent peak_intervals

        track ${group}_nuintervals
        bigDataUrl ${group}_trimmed_peaks.bb
        shortLabel $group non-unique
        longLabel $group: Trimmed peak regions from non-unique reads
        subGroups cline=${ucsc_group[$group]}
        type bigBed 6
        color $color
        parent peak_intervals

intervals_track

    rm -f *.unsorted
    gzip -f *.bed
done
