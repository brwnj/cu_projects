#!/usr/bin/env bash
#BSUB -J ribosome
#BSUB -e ribosome.%J.err
#BSUB -o ribosome.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P kabos_ribosome

<<DOC
ribosome profiling experiment

select[mem>16] rusage[mem=16] span[hosts=1]

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
        cmd="zcat $input_file | fastx_clipper -a $TRIMADAPTER -l $TRIMMINLENGTH -c -n -v -Q33 | fastx_trimmer -Q33 -f 2 | gzip -c > $output_file"
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
        runscript=${sample}_filtering.sh
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
    output_dir=$RESULTS/$sample/alignments
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/accepted_hits.bam
    if [[ ! -f $output_file ]]; then
        cmd="tophat -o $output_dir -p 10 --segment-mismatches $TOPHATSEGMENTMISMATCHES --no-novel-juncs -T --transcriptome-index $TOPHATTRANSCRIPTOMEINDEX $TOPHATREF $input_file"
        bsub -J align -o align.%J.out -e align.%J.err -P $PROJECTID -R "select[mem>16] rusage[mem=16] span[hosts=1]" -n 10 -K $cmd &
    fi
done
wait
