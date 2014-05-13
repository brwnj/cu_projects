#!/usr/bin/env bash
#BSUB -J hitsclip
#BSUB -e hitsclip.%J.err
#BSUB -o hitsclip.%J.out
#BSUB -q normal
#BSUB -n 1
#BSUB -P hits-clip

set -o nounset -o pipefail -o errexit -x
source /vol1/home/brownj/projects/hits-clip/bin/config.sh


# trim the adapter sequence
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=trim
    sample=${SAMPLES[$i]}
    input_file=$DATA/$sample.fastq.gz
    output_dir=$DATA/trimmed
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/${sample}.fastq.gz
    if [[ ! -f $output_file ]]; then
        cmd="python $TRIMSCRIPT -a $TRIMADAPTER -d TCAGTC $input_file | gzip -c > $output_file"
        bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
    fi
done
wait


# align trimmed reads
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    jname=align
	sample=${SAMPLES[$i]}
    input_file=$DATA/trimmed/$sample.fastq.gz
    output_dir=$RESULTS/$sample/alignments/novoalign
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
	output_file_1=$output_dir/alignment_summary.txt
    output_file_2=$output_dir/$sample.bam
    if [[ ! -f $output_file_2 ]]; then
        runscript=${jname}_${sample}.sh
		echo "set -x" > $runscript
        echo "novoalign -d $NOVOIDX -f $input_file -a -o SAM -r A 20 -e 100 -c 10 -n 60 -s 4 -l 16 -k 2> $output_file_1 | samtools view -ShuF4 - | samtools sort -o - $sample.temp -m 8G > $output_file_2" >> $runscript
        echo "samtools index $output_file_2" >> $runscript
        bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -R "select[mem>16] rusage[mem=16] span[hosts=1]" -n 10 -K < $runscript &
    fi
done
wait


# remove duplicates from reads
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    jname=rmdup
	sample=${SAMPLES[$i]}
    input_file=$RESULTS/$sample/alignments/novoalign/$sample.bam
    output_dir=$RESULTS/$sample/alignments/rmdup
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/$sample.bam
    if [[ ! -f $output_file ]]; then
        runscript=${jname}_${sample}.sh
        echo "samtools rmdup -s $input_file $output_file" > $runscript
        echo "samtools index $output_file" >> $runscript
        bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K < $runscript &
    fi
done
wait


# filter the bams
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=filterbams
    sample=${SAMPLES[$i]}

    input_file=$RESULTS/$sample/alignments/novoalign/$sample.bam
    output_dir=$RESULTS/$sample/alignments/novoalign/filtered
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/${sample}.bam
    if [[ ! -f $output_file ]]; then
        cmd="python $FILTERSCRIPT -m all -d GTGTCA -d GTGCCA -d GTGTCT -d GTCTCA -d GTGACA -b $FILTERDOWNSTREAMBASES $input_file $output_file $FASTA"
        bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
    fi

	# also submit a job for the rmdup bam
    input_file=$RESULTS/$sample/alignments/rmdup/$sample.bam
    output_dir=$RESULTS/$sample/alignments/rmdup/filtered
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/${sample}.bam
    if [[ ! -f $output_file ]]; then
        cmd="python $FILTERSCRIPT -m all -d GTGTCA -d GTGCCA -d GTGTCT -d GTCTCA -d GTGACA -b $FILTERDOWNSTREAMBASES $input_file $output_file $FASTA"
        bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
    fi
done
wait


# make stranded bams
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=strandbams
    sample=${SAMPLES[$i]}

    input_file=$RESULTS/$sample/alignments/novoalign/filtered/$sample.bam
    output_dir=$RESULTS/$sample/alignments/novoalign/filtered
	for strand in pos neg; do
	    output_file=$output_dir/${sample}_${strand}.bam
	    if [[ ! -f $output_file ]]; then
	        cmd="samtools view -hb -F 0x10 $input_file > $output_file"
			if [[ "$strand" = *neg* ]]; then
	        	cmd="samtools view -hb -f 0x10 $input_file > $output_file"
	        fi
	        bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
	    fi
	done

	# rmdup
    input_file=$RESULTS/$sample/alignments/rmdup/filtered/$sample.bam
    output_dir=$RESULTS/$sample/alignments/rmdup/filtered
	for strand in pos neg; do
	    output_file=$output_dir/${sample}_${strand}.bam
	    if [[ ! -f $output_file ]]; then
	        cmd="samtools view -hb -F 0x10 $input_file > $output_file"
			if [[ "$strand" = *neg* ]]; then
	        	cmd="samtools view -hb -f 0x10 $input_file > $output_file"
	        fi
	        bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
	    fi
	done
done
wait


# make bedgraphs and bigwigs
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    jname=bedgraph
	sample=${SAMPLES[$i]}

    output_dir_1=$RESULTS/$sample/intervals/bedgraph
	output_dir_2=$RESULTS/$sample/intervals/bigwig
    if [[ ! -d $output_dir_1 ]]; then
        mkdir -p $output_dir_1
    fi
    if [[ ! -d $output_dir_2 ]]; then
        mkdir -p $output_dir_2
    fi

	# stranded interval files
	for strand in pos neg; do
		input_file=$RESULTS/$sample/alignments/novoalign/filtered/${sample}_${strand}.bam
	    symbol="+"
		if [[ "$strand" = *neg* ]]; then
			symbol="-"
		fi
		output_file_1=$output_dir_1/${sample}_${strand}.bedgraph.gz
		output_file_2=$output_dir_2/${sample}_${strand}.bigwig
	    if [[ ! -f $output_file_2 ]]; then
	        runscript=${jname}_${sample}_${strand}.sh
	        echo "bedtools genomecov -strand $symbol -bg -ibam $input_file | bedtools sort -i - > ${output_file_1/.gz}" > $runscript
	        echo "bedGraphToBigWig ${output_file_1/.gz} $SIZES $output_file_2" >> $runscript
			echo "gzip -f ${output_file_1/.gz}" >> $runscript
	        bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K < $runscript &
	    fi
	done

	# run also for rmdup bam
    output_dir_1=$RESULTS/$sample/intervals/bedgraph/rmdup
	output_dir_2=$RESULTS/$sample/intervals/bigwig/rmdup
    if [[ ! -d $output_dir_1 ]]; then
        mkdir -p $output_dir_1
    fi
    if [[ ! -d $output_dir_2 ]]; then
        mkdir -p $output_dir_2
    fi

	# stranded interval files
	for strand in pos neg; do
	    input_file=$RESULTS/$sample/alignments/rmdup/filtered/${sample}_${strand}.bam
	    symbol="+"
		if [[ "$strand" = *neg*  ]]; then
			symbol="-"
		fi
		output_file_1=$output_dir_1/${sample}_${strand}.bedgraph.gz
		output_file_2=$output_dir_2/${sample}_${strand}.bigwig
	    if [[ ! -f $output_file_2 ]]; then
	        runscript=${jname}_${sample}_${strand}.sh
	        echo "bedtools genomecov -strand $symbol -bg -ibam $input_file | bedtools sort -i - > ${output_file_1/.gz}" > $runscript
	        echo "bedGraphToBigWig ${output_file_1/.gz} $SIZES $output_file_2" >> $runscript
			echo "gzip -f ${output_file_1/.gz}" >> $runscript
	        bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K < $runscript &
	    fi
	done
done


# call peaks
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    jname=peaks
	sample=${SAMPLES[$i]}
	for strand in pos neg; do
		input_file=$RESULTS/$sample/alignments/novoalign/filtered/${sample}_${strand}.bam
		output_dir=$RESULTS/$sample/peaks/$strand
		if [[ ! -d $output_dir ]]; then
			mkdir -p $output_dir
		fi
		output_file=$output_dir/${sample}_${strand}_peaks.narrowPeak.gz
		if [[ ! -f $output_file ]]; then
			runscript=${jname}_${sample}_${strand}.sh
			echo "macs2 callpeak -t $input_file --outdir $output_dir -g hs -n ${sample}_${strand} --nomodel --extsize 20 -q $PEAKSNONUNIQUEQ --keep-dup all" > $runscript
			echo "gzip -f ${output_file/.gz}" >> $runscript
			echo "rm -f $output_dir/{*.xls,*.bed}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K < $runscript &
		fi
	done

	# rmdup
	for strand in pos neg; do
		input_file=$RESULTS/$sample/alignments/rmdup/filtered/${sample}_${strand}.bam
		output_dir=$RESULTS/$sample/peaks/$strand
		if [[ ! -d $output_dir ]]; then
			mkdir -p $output_dir
		fi
		output_file=$output_dir/${sample}_${strand}_peaks.narrowPeak.gz
		if [[ ! -f $output_file ]]; then
			runscript=${jname}_${sample}_${strand}.sh
			echo "macs2 callpeak -t $input_file --outdir $output_dir -g hs -n ${sample}_${strand}_rmdup --nomodel --extsize 20 -q $PEAKSUNIQUEQ --keep-dup all" > $runscript
			echo "gzip -f ${output_file/.gz}" >> $runscript
			echo "rm -f $output_dir/{*.xls,*.bed}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K < $runscript &
		fi
	done
done
wait


# build genomedata archive;
tracks=""
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	sample=${SAMPLES[$i]}
	tracks="$tracks -t ${sample}_pos=$RESULTS/$sample/intervals/bedgraph/${sample}_pos.bedgraph.gz"
	tracks="$tracks -t ${sample}_neg=$RESULTS/$sample/intervals/bedgraph/${sample}_neg.bedgraph.gz"
	tracks="$tracks -t ${sample}_pos_rmdup=$RESULTS/$sample/intervals/bedgraph/rmdup/${sample}_pos.bedgraph.gz"
	tracks="$tracks -t ${sample}_neg_rmdup=$RESULTS/$sample/intervals/bedgraph/rmdup/${sample}_neg.bedgraph.gz"
done
if [[ ! -d $GENOMEDATA ]]; then
	for chrom in `ls $FASTAS | sed -rn 's/(chr[0-9XYM]+).*/\1/p'`; do
		jname=genomedata
		runscript=${jname}_${chrom}.sh
		echo "genomedata-load -v --directory-mode -s $FASTAS/$chrom.fa.gz $tracks $GENOMEDATA" > $runscript
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K < $runscript &
	done
else
	echo "Genomedata exists. Load new tracks individually and comment this section out."
	# exit 1
fi

# to add to an existing archive...
# you can get existing track names via `genomedata-info tracknames_continuous $GENOMEDATA`
# genomedata-open-data $GDARCHIVE <new track names>
# zcat <bedgraph.gz> | genomedata-load-data $GDARCHIVE <single track name>
# zcat <bedgraph.gz> | genomedata-load-data $GDARCHIVE <single track name>
# zcat <bedgraph.gz> | genomedata-load-data $GDARCHIVE <single track name>
# genomedata-close-data $GDARCHIVE

wait


# merge peaks and trim across replicates
# for sample in ${replicates[$group]}; do
# for k in "${!REPLICATES[@]}"
# do
#   echo "key  : $k"
#   echo "value: ${REPLICATES[$k]}"
# done
#
# for strand in pos neg; do
#     bedfiles=""
# 	tracks=""
#     output_file_1=$results/${group}_${strand}_peaks.bed.gz
# 	output_file_2=trimmed
#     if [[ ! -f $output_file_1 ]]; then
#         for sample in $replicates; do
#             bedfiles="$bedfiles $RESULTS/$sample/$sample.rmd_${strand}.peaks.qv.passed_filter.bed.gz"
# 			# tracks="$tracks -t ${sample}_filtered.rmd_${strand}"
#         done
#         if [[ $replicate_count == 1 ]]; then
#             gunzip $bedfiles
#             bedClip ${bedfiles/.gz} $SIZES ${out/.gz}
#             gzip -f ${out/.gz} ${bedfiles/.gz}
#         else
#             peaktools-combine-replicates --verbose $bedfiles | bedClip stdin $SIZES ${combined/.gz}
#             gzip -f ${combined/.gz}
#
#         fi
#     fi
#
# 	if [[ ! -f $output_file_2 ]]; then
#         cmd="python ~/projects/hits-clip/bin/scripts/trim_peaks.py -v $tracks $combined $GENOMEDATA | bedtools sort -i - | gzip -c > $trimmed"
#         bsub -J trim -o trim.%J.out -e trim.%J.err -P hits-clip -K $cmd &
# 	fi
#
# 	# then do the same for non rmdup
#
# done
#
# wait


# make a hub with coverage and peaks for samples only
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
    echo "hub $HUBNAME" > $hub
    echo "shortLabel $HUBNAME" >> $hub
    echo "longLabel $HUBNAME" >> $hub
    echo "genomesFile genomes.txt" >> $hub
    echo "email brwnjm@gmail.com" >> $hub
fi

# output for coverage
trackdb=$HUB/$GENOME/trackDb.txt
cat <<coverage_track >$trackdb
track ${HUBNAME}_coverage
compositeTrack on
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
    posbw=$RESULTS/$sample/intervals/bigwig/${sample}_pos.bigwig
    negbw=$RESULTS/$sample/intervals/bigwig/${sample}_neg.bigwig
    posuniqbw=$RESULTS/$sample/intervals/bigwig/rmdup/${sample}_pos.bigwig
    neguniqbw=$RESULTS/$sample/intervals/bigwig/rmdup/${sample}_neg.bigwig

	cp $posbw $negbw $HUB/$GENOME
	cp $posuniqbw $HUB/$GENOME/${sample}_rmdup_pos.bigwig
	cp $neguniqbw $HUB/$GENOME/${sample}_rmdup_neg.bigwig
	# yuck.
	posuniqbw=${sample}_rmdup_pos.bigwig
	neguniqbw=${sample}_rmdup_neg.bigwig

    color=${COLORS[$sample]}
    cat <<coverage_track >>$trackdb
    track $(basename $posbw .bigwig)
    bigDataUrl $(basename $posbw)
    shortLabel $sample coverage POS
    longLabel $sample coverage positive (+) strand
    type bigWig
    parent ${HUBNAME}_coverage
    color $color

    track $(basename $negbw .bigwig)
    bigDataUrl $(basename $negbw)
    shortLabel $sample coverage NEG
    longLabel $sample coverage negative (-) strand
    type bigWig
    parent ${HUBNAME}_coverage
    color $color

    track $(basename $posuniqbw .bigwig)
    bigDataUrl $(basename $posuniqbw)
    shortLabel $sample unique coverage POS
    longLabel $sample unique coverage positive (+) strand
    type bigWig
    parent ${HUBNAME}_coverage
    color $color

    track $(basename $neguniqbw .bigwig)
    bigDataUrl $(basename $neguniqbw)
    shortLabel $sample unique coverage NEG
    longLabel $sample unique coverage negative (-) strand
    type bigWig
    parent ${HUBNAME}_coverage
    color $color

coverage_track
done

# intervals
cat <<intervals_track >>$trackdb
track peak_intervals
compositeTrack on
shortLabel Peak Intervals
longLabel Peak Intervals
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
    color=${COLORS[$sample]}

    # non-unique
	bigbed=$HUB/$GENOME/${sample}_peaks.narrowPeak.bb
	if [[ ! -f $bigbed ]]; then

	    np_us=${sample}.narrowPeak.unsorted
	    np=${sample}.narrowPeak
	    bb=${np/.narrowPeak/.bb}

		negnp=$RESULTS/$sample/peaks/neg/${sample}_neg_peaks.narrowPeak.gz
		posnp=$RESULTS/$sample/peaks/pos/${sample}_pos_peaks.narrowPeak.gz

	    awk -t -cbed '{if($5>1000){$5=1000}; $6="-"; print}' $negnp > $np_us
	    awk -t -cbed '{if($5>1000){$5=1000}; $6="+"; print}' $posnp >> $np_us
	    bedtools sort -i $np_us > $np
	    bedToBigBed -type=bed6+4 -as=$fields $np $SIZES $bigbed
	    rm -f $np_us $np

	    cat <<intervals_track >>$trackdb
        track ${bb/.bb}
        bigDataUrl $bb
        shortLabel ${bb/.bb} non-unique
        longLabel $sample: non-unique peaks
        type bigBed 6 +
        color $color
        parent peak_intervals

intervals_track

	fi

    # unique
	bigbeg=$HUB/$GENOME/${sample}_rmdup_peaks.narrowPeak.bb
	if [[ ! -f $bigbed ]]; then

	    np_us=${sample}_rmdup.narrowPeak.unsorted
	    np=${sample}_rmdup.narrowPeak
	    bb=${np/.narrowPeak/.bb}

		negnp=$RESULTS/$sample/peaks/neg/${sample}_neg_rmdup_peaks.narrowPeak.gz
		posnp=$RESULTS/$sample/peaks/pos/${sample}_pos_rmdup_peaks.narrowPeak.gz

	    awk -t -cbed '{if($5>1000){$5=1000}; $6="-"; print}' $negnp > $np_us
	    awk -t -cbed '{if($5>1000){$5=1000}; $6="+"; print}' $posnp >> $np_us
	    bedtools sort -i $np_us > $np
	    bedToBigBed -type=bed6+4 -as=$fields $np $SIZES $bigbed
	    rm -f $np_us $np

	    cat <<intervals_track >>$trackdb
	        track ${bb/.bb}
	        bigDataUrl $bb
	        shortLabel ${bb/.bb} unique
	        longLabel $sample: unique peaks
	        type bigBed 6 +
	        color $color
	        parent peak_intervals

intervals_track

	fi

done
