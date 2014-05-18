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
jname=genomedata
if [[ ! -d $GENOMEDATA ]]; then
	for chrom in `ls $FASTAS | sed -rn 's/(chr[0-9XYM]+).*/\1/p'`; do
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
src=/vol1/home/brownj/projects/hits-clip/results/20140317/peaktools/peaktools
jname=combineandtrim
for group in "${!REPLICATES[@]}"; do
	samples=${REPLICATES[$group]}
	for strand in pos neg; do

		bedfiles=""
		tracks=""

		outdir=$RESULTS/$group/peaks/$strand/postprocessing
		if [[ ! -d $outdir ]]; then
			mkdir -p $outdir
		fi

		output_file_1=$outdir/${group}_${strand}_combined.bed.gz
		output_file_2=$outdir/${group}_${strand}_trimmed.bed.gz

		if [[ ! -f $output_file_1 && ! -f $output_file_2 ]]; then
			for sample in $samples; do
				bedfiles="$bedfiles $RESULTS/$sample/peaks/$strand/${sample}_${strand}_peaks.narrowPeak.gz"
				tracks="$tracks -t ${sample}_${strand}"
			done

			runscript=${jname}_${sample}_${strand}.sh
			echo "python $src/combine_replicates.py --verbose $bedfiles | bedClip stdin $SIZES ${output_file_1/.gz}" > $runscript
			echo "gzip -f ${output_file_1/.gz}" >> $runscript
			echo "python $SRC/trim_peaks.py -v $tracks $output_file_1 $GENOMEDATA | bedtools sort -i - | gzip -c > $output_file_2" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K < $runscript &
		fi

		# then do the same for rmdup
		bedfiles=""
		tracks=""

		if [[ ! -d $outdir ]]; then
			mkdir -p $outdir
		fi

		output_file_1=$outdir/${group}_${strand}_rmdup_combined.bed.gz
		output_file_2=$outdir/${group}_${strand}_rmdup_trimmed.bed.gz

		if [[ ! -f $output_file_1 && ! -f $output_file_2 ]]; then
			for sample in $samples; do
				bedfiles="$bedfiles $RESULTS/$sample/peaks/$strand/${sample}_${strand}_rmdup_peaks.narrowPeak.gz"
				tracks="$tracks -t ${sample}_${strand}_rmdup"
			done

			runscript=${jname}_${sample}_${strand}_rmdup.sh
			echo "python $src/combine_replicates.py --verbose $bedfiles | bedClip stdin $SIZES ${output_file_1/.gz}" > $runscript
			echo "gzip -f ${output_file_1/.gz}" >> $runscript
			echo "python $SRC/trim_peaks.py -v $tracks $output_file_1 $GENOMEDATA | bedtools sort -i - | gzip -c > $output_file_2" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K < $runscript &
		fi

	done
done


# trim the peaks from singletons
jname=trim
for (( i = 0; i < ${#SINGLETONS[@]}; i++ )); do
	sample=${SINGLETONS[$i]}

	for strand in pos neg; do

		track="-t ${sample}_${strand}"
		input_file=$RESULTS/$sample/peaks/$strand/${sample}_${strand}_peaks.narrowPeak.gz
		outdir=$RESULTS/$sample/peaks/$strand/postprocessing

		if [[ ! -d $outdir ]]; then
			mkdir -p $outdir
		fi

		output_file=$outdir/${sample}_${strand}_trimmed.bed.gz

		if [[ ! -f $output_file ]]; then
			cmd="python $SRC/trim_peaks.py -v $track $input_file $GENOMEDATA | bedtools sort -i - | gzip -c > $output_file"
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
		fi

		# then do the same for rmdup
		track="-t ${sample}_${strand}_rmdup"
		input_file=$RESULTS/$sample/peaks/$strand/${sample}_${strand}_rmdup_peaks.narrowPeak.gz

		output_file=$outdir/${sample}_${strand}_rmdup_trimmed.bed.gz

		if [[ ! -f $output_file ]]; then
			cmd="python $SRC/trim_peaks.py -v $track $input_file $GENOMEDATA | bedtools sort -i - | gzip -c > $output_file"
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
		fi

	done
done
wait


# map seeds

#     POSPEAK="$DATA/$SAMPLE.pos.peaks.trimmed.bed.gz"
#     NEGPEAK="$DATA/$SAMPLE.neg.peaks.trimmed.bed.gz"
#
#     SIF="$SAMPLE.network.unfiltered.sif"
#     FUNCTION="$SAMPLE.function.noa"
#
#     RUN=$SAMPLE.mirnaseeds.sh
#
#
# for STRAND in pos neg; do
#     for (( i = 0; i < \${#FASTAS[@]}; i++ )); do
#         SEED_FASTA=\${FASTAS[\$i]}
#         SEED_LENGTH=\${LENGTHS[\$i]}
#
#         POSNOTSIF=\"$SAMPLE.pos.mirna.\$SEED_LENGTH.notSIF.gz\"
#         NEGNOTSIF=\"$SAMPLE.neg.mirna.\$SEED_LENGTH.notSIF.gz\"
#
#         if [ \"\$STRAND\" == \"pos\" ]; then
#             python $SRC/aggregate_peaks.py $POSPEAK $ANNOTATION \$SEED_FASTA $GDARCHIVE | gzip -c > \$POSNOTSIF
#             awk -v len=\"\$SEED_LENGTH\" -c header '{split(\$1,seed,\"|\"); print seed[1]\"\t\"len\"\t\"\$5}' \$POSNOTSIF > $SIF
#         else
#             python $SRC/aggregate_peaks.py $NEGPEAK $ANNOTATION \$SEED_FASTA $GDARCHIVE | gzip -c > \$NEGNOTSIF
#             awk -v len=\"\$SEED_LENGTH\" -c header '{split(\$1,seed,\"|\"); print seed[1]\"\t\"len\"\t\"\$5}' \$NEGNOTSIF >> $SIF
#         fi
#
#     done
# done
#
# echo "functionalCategory" > $FUNCTION
# awk -c header '{print \$1\" = miRNA\n\"\$3\" = Gene\"}' $SIF | sort | uniq >> $FUNCTION
#
# done


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
	bigbed=$HUB/$GENOME/${sample}_peaks.bb
	if [[ ! -f $bigbed ]]; then

	    np_us=${sample}_peaks.narrowPeak.unsorted
	    np=${sample}_peaks.narrowPeak
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
	bigbed=$HUB/$GENOME/${sample}_rmdup_peaks.bb
	if [[ ! -f $bigbed ]]; then

	    np_us=${sample}_rmdup_peaks.narrowPeak.unsorted
	    np=${sample}_rmdup_peaks.narrowPeak
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
