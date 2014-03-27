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

# sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

fasta=$HOME/ref/hg19/hg19.fa
fasta_dir=$HOME/ref/hg19/fa_per_chr
sizes=$HOME/ref/hg19/hg19.sizes
results=/vol1/home/brownj/projects/hits-clip/results/20140317

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
# for track in `genomedata-info tracknames_continuous $GENOMEDATA`; do
#
#     symbol="+"
#     if [[ "${track#*pos}" == "$track" ]]; then
#         symbol="-"
#     fi
#
#     echo "/vol1/software/modules-python/python/2.7.2/bin/peaktools-identify-peaks -t $track -w 30 -v -s $symbol $GENOMEDATA | bedtools sort -i - | gzip -c > $track.peaks.bed.gz" | bsez -j peaks -p hits-clip
#     echo "/vol1/software/modules-python/python/2.7.2/bin/peaktools-identify-peaks -t $track -w 30 -v -s $symbol --shuffle-data $GENOMEDATA | bedtools sort -i - | gzip -c > $track.shuffle.peaks.bed.gz" | bsez -j null_peaks -p hits-clip
# done

# qvalues for peaks
# cutoff=.02
# for track in `genomedata-info tracknames_continuous $GENOMEDATA`; do
#     sample=$(echo ${track/_*} | cut -f1 -d.)
#     real=$track.peaks.bed.gz
#     null=$track.shuffle.peaks.bed.gz
#     qvalues=$track.peaks.unfiltered_qvalues.bed.gz
#     echo "/vol1/software/modules-python/python/2.7.2/bin/peaktools-qvalues -v $real $null | gzip -c > $qvalues" | bsez -j qvalues -p hits-clip
#     # awk -v CO=${cutoff} '\$7<CO'
# done


#################
# PEAKS
#################
# rm *.neg.bam
# rm *.pos.bam
# for bam in *.bam; do
#     # split the bam into positive and negative strands
#     negbam=${bam/.bam/.neg.bam}
#     if [[ ! -f $negbam ]]; then
#         samtools view -hb -f 0x10 $bam > $negbam
#     fi
#     # call peaks on the stranded bams
#     if [[ ! -f ${negbam/.bam}_summits.bed ]]; then
#         macs2 callpeak -t $negbam -g hs -n ${negbam/.bam} --nomodel -s 28 --extsize 15 -q 0.01
#     fi
#     # extend peak regions to allow more likely overlap
#     bedtools slop -i ${negbam/.bam}_summits.bed -g $sizes -b 50 | gzip -c > ${negbam/.bam/.extended_peaks.bed.gz}
#
#     # positive strand
#     posbam=${bam/.bam/.pos.bam}
#     if [[ ! -f $posbam ]]; then
#         samtools view -hb -F 0x10 $bam > $posbam
#     fi
#     if [[ ! -f ${posbam/.bam}_summits.bed ]]; then
#         macs2 callpeak -t $posbam -g hs -n ${posbam/.bam} --nomodel -s 28 --extsize 15 -q 0.01
#     fi
#     bedtools slop -i ${posbam/.bam}_summits.bed -g $sizes -b 50 | gzip -c > ${posbam/.bam/.extended_peaks.bed.gz}
# done


# merge peaks across replicates
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
ucsc_group=([PK11]="MCF7"
            [PK31]="MCF7"
            [PK51]="MCF7"
            [PK21]="BT474"
            [PK41]="BT474"
            [PK52]="BT474"
            [PK24]="MDA231"
            [PK42]="MDA231"
            [PK54]="MDA231")

# for group in MCF7 BT474 MDA231; do
#     for strand in pos neg; do
#         bedfiles=""
#         out=${group}_${strand}_rmd_peaks.bed.gz
#         for sample in ${replicates[$group]}; do
#             bedfiles="$bedfiles ${sample}_filtered.rmd.${strand}.extended_peaks.bed.gz"
#         done
#         toclip=${out/.gz/.clipme}
#         peaktools-combine-replicates --verbose $bedfiles > $toclip
#         bedClip $toclip $sizes ${out/.gz}
#         gzip -f ${out/.gz}
#         rm -f $toclip
#
#         # also run bams containing duplicates
#         bedfiles=""
#         out=${group}_${strand}_peaks.bed.gz
#         for sample in ${replicates[$group]}; do
#             bedfiles="$bedfiles ${sample}_filtered.${strand}.extended_peaks.bed.gz"
#         done
#         toclip=${out/.gz/.clipme}
#         peaktools-combine-replicates --verbose $bedfiles > $toclip
#         bedClip $toclip $sizes ${out/.gz}
#         gzip -f ${out/.gz}
#         rm -f $toclip
#     done
# done

# trim the peaks
for group in MCF7 BT474 MDA231; do
    for strand in pos neg; do
        tracks=""
        out=${group}_${strand}_rmd_peaks.bed.gz
        untrimmed=$results/${group}_${strand}_rmd_peaks.bed.gz
        trimmed=$results/${group}_${strand}_rmd_trimmed_peaks.bed.gz

        for sample in ${replicates[$group]}; do
            tracks="$tracks -t ${sample}_filtered.rmd_${strand}"
        done
        echo "python ~/projects/hits-clip/bin/scripts/trim_peaks.py -v $tracks $untrimmed $GENOMEDATA | bedtools sort -i - | gzip -c > $trimmed" | bsez -j trim -p hits-clip

        #
        tracks=""
        out=${group}_${strand}_peaks.bed.gz
        untrimmed=$results/${group}_${strand}_peaks.bed.gz
        trimmed=$results/${group}_${strand}_trimmed_peaks.bed.gz

        for sample in ${replicates[$group]}; do
            tracks="$tracks -t ${sample}_filtered_${strand}"
        done
        echo "python ~/projects/hits-clip/bin/scripts/trim_peaks.py -v $tracks $untrimmed $GENOMEDATA | bedtools sort -i - | gzip -c > $trimmed" | bsez -j trim -p hits-clip

    done
done


# add to some hub somewhere
# coverage, peaks.

# HUB=hub
# GENOME=hg19
#
# if [[ ! -d $HUB/$GENOME ]]; then
#     mkdir -p $HUB/$GENOME
# fi
#
# # genomes.txt
# if [[ ! -f $HUB/genomes.txt ]]; then
#     genomes=$HUB/genomes.txt
#     echo "genome $GENOME" > $genomes
#     echo "trackDb $GENOME/trackDb.txt" >> $genomes
# fi
#
# # hub.txt
# if [[ ! -f $HUB/hub.txt ]]; then
#     hub=$HUB/hub.txt
#     echo "hub re_hitsclip" > $hub
#     echo "shortLabel re-HITS-CLIP" >> $hub
#     echo "longLabel re-HITS-CLIP" >> $hub
#     echo "genomesFile genomes.txt" >> $hub
#     echo "email brwnjm@gmail.com" >> $hub
# fi
#
# trackdb=$HUB/$GENOME/trackDb.txt
#
# # output for coverage
# cat <<coverage_track >$trackdb
# track re_hitsclip_coverage
# compositeTrack on
# subGroup1 cline CellLine MDA231=MDA231 BT474=BT474 MCF7=MCF7
# shortLabel Coverage
# longLabel Coverage
# maxHeightPixels 50:20:15
# type bigWig
# configurable on
# autoScale on
#
# coverage_track
#
# for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
#     sample=${SAMPLES[$i]}
#
#     posbg=${sample}_filtered_pos.bedgraph.gz
#     negbg=${sample}_filtered_neg.bedgraph.gz
#     posuniqbg=${sample}_filtered.rmd_pos.bedgraph.gz
#     neguniqbg=${sample}_filtered.rmd_neg.bedgraph.gz
#
#     gunzip $posbg $negbg $posuniqbg $neguniqbg
#
#     posbw=${sample}_filtered_pos.bw
#     negbw=${sample}_filtered_neg.bw
#     posuniqbw=${sample}_filtered.rmd_pos.bw
#     neguniqbw=${sample}_filtered.rmd_neg.bw
#
#     bedGraphToBigWig ${posbg/.gz} $sizes $posbw
#     bedGraphToBigWig ${negbg/.gz} $sizes $negbw
#     bedGraphToBigWig ${posuniqbg/.gz} $sizes $posuniqbw
#     bedGraphToBigWig ${neguniqbg/.gz} $sizes $neguniqbw
#
#     gzip ${posbg/.gz} ${negbg/.gz} ${posuniqbg/.gz} ${neguniqbg/.gz}
#
#     mv $posbw $negbw $posuniqbw $neguniqbw $HUB/$GENOME
#
#     color=${colors[$sample]}
#     cat <<coverage_track >>$trackdb
#     track ${posbw/.bw}
#     bigDataUrl $posbw
#     subGroups cline=${ucsc_group[$sample]}
#     shortLabel $sample coverage POS
#     longLabel $sample coverage positive (+) strand
#     type bigWig
#     parent re_hitsclip_coverage
#     color $color
#
#     track ${negbw/.bw}
#     bigDataUrl $negbw
#     subGroups cline=${ucsc_group[$sample]}
#     shortLabel $sample coverage NEG
#     longLabel $sample coverage negative (-) strand
#     type bigWig
#     parent re_hitsclip_coverage
#     color $color
#
#     track ${posuniqbw/.bw}
#     bigDataUrl $posuniqbw
#     subGroups cline=${ucsc_group[$sample]}
#     shortLabel $sample unique coverage POS
#     longLabel $sample unique coverage positive (+) strand
#     type bigWig
#     parent re_hitsclip_coverage
#     color $color
#
#     track ${neguniqbw/.bw}
#     bigDataUrl $neguniqbw
#     subGroups cline=${ucsc_group[$sample]}
#     shortLabel $sample unique coverage NEG
#     longLabel $sample unique coverage negative (-) strand
#     type bigWig
#     parent re_hitsclip_coverage
#     color $color
#
# coverage_track
# done
#
# # intervals
# cat <<intervals_track >>$trackdb
# track peak_intervals
# compositeTrack on
# shortLabel Peak Intervals
# longLabel Peak Intervals
# type bigBed 6
#
# intervals_track
#
# awk -t -cbed '{print $1, $2, $3, "BT474_rmd", "0", "+"}' BT474_pos_rmd_trimmed_peaks.bed.gz >> BT474_rmd_trimmed_peaks.unsorted
# awk -t -cbed '{print $1, $2, $3, "BT474_rmd", "0", "-"}' BT474_neg_rmd_trimmed_peaks.bed.gz >> BT474_rmd_trimmed_peaks.unsorted
# bedtools sort -i BT474_rmd_trimmed_peaks.unsorted | bedtools merge -s -i - | bedtools sort -i - | awk -t '{print $1, $2, $3, "BT474_rmd", 0, $6}' > BT474_rmd_trimmed_peaks.bed
# bedToBigBed -type=bed6 BT474_rmd_trimmed_peaks.bed ~/ref/hg19/hg19.sizes BT474_rmd_trimmed_peaks.bb
# mv BT474_rmd_trimmed_peaks.bb $HUB/$GENOME
# awk -t -cbed '{print $1, $2, $3, "BT474", "0", "+"}' BT474_pos_trimmed_peaks.bed.gz >> BT474_trimmed_peaks.unsorted
# awk -t -cbed '{print $1, $2, $3, "BT474", "0", "-"}' BT474_neg_trimmed_peaks.bed.gz >> BT474_trimmed_peaks.unsorted
# bedtools sort -i BT474_trimmed_peaks.unsorted | bedtools merge -s -i - | bedtools sort -i - | awk -t '{print $1, $2, $3, "BT474", 0, $6}' > BT474_trimmed_peaks.bed
# bedToBigBed -type=bed6 BT474_trimmed_peaks.bed ~/ref/hg19/hg19.sizes BT474_trimmed_peaks.bb
# mv BT474_trimmed_peaks.bb $HUB/$GENOME
# awk -t -cbed '{print $1, $2, $3, "MCF7_rmd", "0", "+"}' MCF7_pos_rmd_trimmed_peaks.bed.gz >> MCF7_rmd_trimmed_peaks.unsorted
# awk -t -cbed '{print $1, $2, $3, "MCF7_rmd", "0", "-"}' MCF7_neg_rmd_trimmed_peaks.bed.gz >> MCF7_rmd_trimmed_peaks.unsorted
# bedtools sort -i MCF7_rmd_trimmed_peaks.unsorted | bedtools merge -s -i - | bedtools sort -i - | awk -t '{print $1, $2, $3, "MCF7_rmd", 0, $6}' > MCF7_rmd_trimmed_peaks.bed
# bedToBigBed -type=bed6 MCF7_rmd_trimmed_peaks.bed ~/ref/hg19/hg19.sizes MCF7_rmd_trimmed_peaks.bb
# mv MCF7_rmd_trimmed_peaks.bb $HUB/$GENOME
# awk -t -cbed '{print $1, $2, $3, "MCF7", "0", "+"}' MCF7_pos_trimmed_peaks.bed.gz >> MCF7_trimmed_peaks.unsorted
# awk -t -cbed '{print $1, $2, $3, "MCF7", "0", "-"}' MCF7_neg_trimmed_peaks.bed.gz >> MCF7_trimmed_peaks.unsorted
# bedtools sort -i MCF7_trimmed_peaks.unsorted | bedtools merge -s -i - | bedtools sort -i - | awk -t '{print $1, $2, $3, "MCF7", 0, $6}' > MCF7_trimmed_peaks.bed
# bedToBigBed -type=bed6 MCF7_trimmed_peaks.bed ~/ref/hg19/hg19.sizes MCF7_trimmed_peaks.bb
# mv MCF7_trimmed_peaks.bb $HUB/$GENOME
# awk -t -cbed '{print $1, $2, $3, "MDA231_rmd", "0", "+"}' MDA231_pos_rmd_trimmed_peaks.bed.gz >> MDA231_rmd_trimmed_peaks.unsorted
# awk -t -cbed '{print $1, $2, $3, "MDA231_rmd", "0", "-"}' MDA231_neg_rmd_trimmed_peaks.bed.gz >> MDA231_rmd_trimmed_peaks.unsorted
# bedtools sort -i MDA231_rmd_trimmed_peaks.unsorted | bedtools merge -s -i - | bedtools sort -i - | awk -t '{print $1, $2, $3, "MDA231_rmd", 0, $6}' > MDA231_rmd_trimmed_peaks.bed
# bedToBigBed -type=bed6 MDA231_rmd_trimmed_peaks.bed ~/ref/hg19/hg19.sizes MDA231_rmd_trimmed_peaks.bb
# mv MDA231_rmd_trimmed_peaks.bb $HUB/$GENOME
# awk -t -cbed '{print $1, $2, $3, "MDA231", "0", "+"}' MDA231_pos_trimmed_peaks.bed.gz >> MDA231_trimmed_peaks.unsorted
# awk -t -cbed '{print $1, $2, $3, "MDA231", "0", "-"}' MDA231_neg_trimmed_peaks.bed.gz >> MDA231_trimmed_peaks.unsorted
# bedtools sort -i MDA231_trimmed_peaks.unsorted | bedtools merge -s -i - | bedtools sort -i - | awk -t '{print $1, $2, $3, "MDA231", 0, $6}' > MDA231_trimmed_peaks.bed
# bedToBigBed -type=bed6 MDA231_trimmed_peaks.bed ~/ref/hg19/hg19.sizes MDA231_trimmed_peaks.bb
# mv MDA231_trimmed_peaks.bb $HUB/$GENOME
# gzip -f *.bed
# rm -f *.unsorted
#
# for group in BT474 MCF7 MDA231; do
#     color=${colors[$group]}
#     cat <<intervals_track >>$trackdb
#         track ${group}_uintervals
#         bigDataUrl ${group}_rmd_trimmed_peaks.bb
#         shortLabel $group unique
#         longLabel $group: Trimmed peak regions from unique reads
#         type bigBed 6
#         color $color
#         parent peak_intervals
#
#         track ${group}_nuintervals
#         bigDataUrl ${group}_trimmed_peaks.bb
#         shortLabel $group non-unique
#         longLabel $group: Trimmed peak regions from non-unique reads
#         type bigBed 6
#         color $color
#         parent peak_intervals
#
# intervals_track
# done
