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
        echo "novoalign -d $NOVOIDX -f $input_file -a -o SAM -r A 20 -e 100 -c 10 -s 4 -l 16 -k 2> $output_file_1 | samtools view -ShuF4 - | samtools sort -o - $sample.temp -m 8G > $output_file_2" > $runscript
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

    input_file=$RESULTS/$sample/alignments/novoalign/rmdup/$sample.bam
    output_dir=$RESULTS/$sample/alignments/novoalign/filtered/
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





exit





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


exit



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
		output_file=$output_dir/${sample}_peaks.narrowPeak.gz
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
		output_file=$output_dir/${sample}_peaks.narrowPeak.gz
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
fi

# to add to an existing archive...
# you can get existing track names via `genomedata-info tracknames_continuous $GENOMEDATA`
# genomedata-open-data $GDARCHIVE <new track names>
# zcat <bedgraph.gz> | genomedata-load-data $GDARCHIVE <single track name>
# zcat <bedgraph.gz> | genomedata-load-data $GDARCHIVE <single track name>
# zcat <bedgraph.gz> | genomedata-load-data $GDARCHIVE <single track name>
# genomedata-close-data $GDARCHIVE

wait


exit

# merge peaks across replicates
for k in "${!REPLICATES[@]}"
do
  echo "key  : $k"
  echo "value: ${array[$k]}"
done

for strand in pos neg; do
    bedfiles=""
    out=$results/$group.$strand.peaks.bed.gz
    if [[ ! -f $out ]]; then
        for sample in $replicates; do
            bedfiles="$bedfiles $RESULTS/$sample/$sample.rmd_${strand}.peaks.qv.passed_filter.bed.gz"
        done
        if [[ $replicate_count == 1 ]]; then
            gunzip $bedfiles
            bedClip ${bedfiles/.gz} $SIZES ${out/.gz}
            gzip -f ${out/.gz}
            gzip ${bedfiles/.gz}
        else
            toclip=${out/.gz/.clipme}
            peaktools-combine-replicates --verbose $bedfiles > $toclip
            bedClip $toclip $SIZES $out
            gzip $out
            rm -f $toclip
        fi
    fi
done


# trim peaks
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
wait