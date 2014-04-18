#!/usr/bin/env bash
#BSUB -J ribosome
#BSUB -e ribosome.%J.err
#BSUB -o ribosome.%J.out
#BSUB -q normal
#BSUB -R "select[mem>1] rusage[mem=1] span[hosts=1]"
#BSUB -n 1
#BSUB -P kabos_ribosome

<<DOC
ribosome profiling experiment
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/kabos_ribosome/bin/config.sh


# trim the sequences using fastx toolkit
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}
    input_file=$DATA/$sample.fastq.gz
    output_dir=$DATA/trimmed
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/${sample}.fastq.gz
    if [[ ! -f $output_file ]]; then
        cmd="zcat $input_file | fastx_clipper -a $TRIMADAPTER -l $TRIMMINLENGTH -c -n -v -Q33 | fastx_trimmer -Q33 -f 1 | gzip -c > $output_file"
        bsub -J trim -o trim.%J.out -e trim.%J.err -P $PROJECTID -K $cmd &
    fi
done
wait


# filter sequences using bowtie and rRNA database
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}
    input_file=$DATA/trimmed/$sample.fastq.gz
    output_dir=$DATA/filtered
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/$sample.fastq.gz
    if [[ ! -f $output_file ]]; then
        tmp_file=$sample.tmp
        runscript=filtering_${sample}.sh
        echo "zcat $input_file | bowtie --un $tmp_file --seedlen $FILTERSEEDLEN -p 10 $FILTERINDEX - > /dev/null" > $runscript
        echo "gzip $tmp_file" >> $runscript
        echo "mv $tmp_file.gz $output_file" >> $runscript
        bsub -J filter -o filter.%J.out -e filter.%J.err -P $PROJECTID -R "select[mem>16] rusage[mem=16] span[hosts=1]" -n 10 -K < $runscript &
    fi
done
wait


# align reads using tophat as per methods paper
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}
    input_file=$DATA/filtered/$sample.fastq.gz
    output_dir=$RESULTS/$sample/alignments/star
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/$sample.bam
    if [[ ! -f $output_file ]]; then
        cpus=12
        tmp_file=$output_dir/${sample}_Aligned.out.sam
        runscript=star_${sample}.sh
        echo "STAR --runThreadN $cpus --genomeDir $STARGENOMEDIR --readFilesIn $input_file --readFilesCommand zcat --outFileNamePrefix $output_dir/${sample}_ --outFilterMultimapNmax 5 --sjdbGTFfile $HG19GTF" > $runscript
        echo "samtools view -ShuF4 $tmp_file | samtools sort -@ $cpus -m 16G - ${output_file/.bam}" >> $runscript
        echo "samtools index $output_file" >> $runscript
        echo "rm $tmp_file" >> $runscript
        bsub -J align -o align.%J.out -e align.%J.err -P $PROJECTID -R "select[mem>16] rusage[mem=16] span[hosts=1]" -n $cpus -K < $runscript &
    fi
done
wait


# create bedgraphs for ribosome footprint QC
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}
    input_file=$RESULTS/$sample/alignments/star/$sample.bam
    output_dir=$RESULTS/$sample/interval_files
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/${sample}_5p_coverage.bedgraph.gz
    if [[ ! -f $output_file ]]; then
        cmd="bedtools genomecov -5 -bg -ibam $input_file | bedtools sort -i - | gzip -c > $output_file"
        bsub -J gcov -o gcov.%J.out -e gcov.%J.err -P $PROJECTID -K $cmd &
    fi
done
wait


# output ribosome footprint stats
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}
    # only test each fraction against totalRNA
    if [[ "$sample" = *fraction* ]]; then
        celltype=${sample%_*}
        input_file_1=$RESULTS/$sample/interval_files/${sample}_5p_coverage.bedgraph.gz
        input_file_2=$RESULTS/${celltype}_totalRNA/interval_files/${celltype}_totalRNA_5p_coverage.bedgraph.gz
        output_dir=$RESULTS/$sample/postprocessing/ribosome_footprint/
        if [[ ! -d $output_dir ]]; then
            mkdir -p $output_dir
        fi
        output_file=$output_dir/$sample.txt
        script=/vol1/home/brownj/projects/kabos_ribosome/bin/ribosome_footprints.py
        if [[ ! -f $output_file ]]; then
            cmd="python $script $input_file_1 $input_file_2 $HG19STARTCODONS > $output_file"
            bsub -J footprint -o fp.%J.out -e fp.%J.err -P $PROJECTID -K $cmd &
        fi
    fi
done
wait


# make a hub; coverage tracks
if [[ ! -d $HUB/$GENOME ]]; then
    mkdir -p $HUB/$GENOME
fi

# genomes.txt
if [[ ! -f $GENOMES ]]; then
    echo "genome $GENOME" > $GENOMES
    echo "trackDb $GENOME/trackDb.txt" >> $GENOMES
fi

# hub.txt
if [[ ! -f $HUB/hub.txt ]]; then
    hub=$HUB/hub.txt
    echo "hub $HUBNAME" > $hub
    echo "shortLabel $HUBLABEL" >> $hub
    echo "longLabel $HUBLABEL" >> $hub
    echo "genomesFile genomes.txt" >> $hub
    echo "email $HUBEMAIL" >> $hub
fi

# make bigwigs
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}
    input_file=$RESULTS/$sample/alignments/star/$sample.bam
    output_dir=$HUB/$GENOME
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    for strand in pos neg; do
        output_file=$output_dir/${sample}_$strand.bw
        if [[ ! -f $output_file ]]; then
            symbol="-"
            if [[ $strand == *pos* ]]; then
                symbol="+"
            fi
            tmp=${sample}_$strand.tmp
            runscript=bam2bw_${sample}_${strand}.sh
            echo "bedtools genomecov -strand $symbol -bg -ibam $input_file | bedtools sort -i - > $tmp" > $runscript
            echo "bedGraphToBigWig $tmp $SIZES $output_file" >> $runscript
            echo "rm $tmp" >> $runscript
            bsub -J bam2bw -o bam2bw.%J.out -e bam2bw.%J.err -P $PROJECTID -K < $runscript &
        fi
    done
done
wait


# output for coverage
cat <<coverage_track >$TRACKDB
track coverage
compositeTrack on
subGroup1 sgroup SampleGroup MCF7=MCF7 PK12=PK12
shortLabel Coverage
longLabel Coverage
maxHeightPixels 50:20:15
type bigWig
configurable on
autoScale on

coverage_track

for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}

    posbw=${sample}_pos.bw
    negbw=${sample}_neg.bw

    color=${COLORS[$sample]}
    cat <<coverage_track >>$TRACKDB
    track ${posbw/.bw}
    bigDataUrl $posbw
    subGroups sgroup=${SUBGROUP1[$sample]}
    shortLabel $sample coverage POS
    longLabel $sample coverage positive (+) strand
    type bigWig
    parent coverage
    color $color

    track ${negbw/.bw}
    bigDataUrl $negbw
    subGroups sgroup=${SUBGROUP1[$sample]}
    shortLabel $sample coverage NEG
    longLabel $sample coverage negative (-) strand
    type bigWig
    parent coverage
    color $color

coverage_track
done


# get counts
input_file=$RESULTS/*/alignments/star/*.bam
output_dir=$POSTPROCESSING/feature_counts
if [[ ! -d $output_dir ]]; then
    mkdir -p $output_dir
fi
output_file=$output_dir/counts.txt
if [[ ! -f $output_file ]]; then
    cmd="featureCounts -a $HG19GTF -o $output_file -T 12 $input_file"
    bsub -J counts -o counts.%J.out -e counts.%J.err -P $PROJECTID -R "span[hosts=1]" -n 12 -K $cmd &
fi
wait

# do the comparisons
