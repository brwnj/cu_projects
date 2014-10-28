#!/usr/bin/env bash
#BSUB -J ichiro
#BSUB -e ichiro.%J.%I.err
#BSUB -o ichiro.%J.%I.out
#BSUB -q normal
#BSUB -R "span[hosts=1]"
#BSUB -P koch

set -o nounset -o pipefail -o errexit -x

samples=(
96B1_AGTTCC_L007
96B2_ATGTCA_L007
96B3_TGACCA_L008
LB6_1_CAGATC_L008
LUDB_8_9_GTCCGC_L007
LUDB_ACAGTG_L008
LUL_GCCAAT_L008
RB8_9_CCGTCC_L007
)

results=$HOME/projects/ichiro_rnaseq/results/common
data=$HOME/projects/ichiro_rnaseq/data/common
pid=ichiro_rnaseq
stargenomedir=$HOME/ref/hg19
hg19gtf=$HOME/ref/hg19/hg19.gtf
hg19reference=$HOME/ref/hg19/hg19.fa
snpeff=$HOME/opt/snpeff/snpEff.jar
snpeffconfig=$HOME/opt/snpeff/snpEff.config


# quality trim the reads
for (( i = 0; i < ${#samples[@]}; i++ )); do
    sample=${samples[$i]}
    input_file_1=$data/${sample}_R1_001.fastq.gz
    input_file_2=$data/${sample}_R2_001.fastq.gz
    output_dir=$data/trimmed
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file_1=$output_dir/${sample}_R1_001.fastq.gz
    if [[ ! -f $output_file_1 ]]; then
        cmd="seqtk trimfq $input_file_1 | gzip -c > $output_file_1"
        bsub -J trim -o trim.%J.out -e trim.%J.err -P $pid -K $cmd &
    fi
    output_file_2=$output_dir/${sample}_R2_001.fastq.gz
    if [[ ! -f $output_file_2 ]]; then
        cmd="seqtk trimfq $input_file_2 | gzip -c > $output_file_2"
        bsub -J trim -o trim.%J.out -e trim.%J.err -P $pid -K $cmd &
    fi
done
wait


# align reads using tophat as per methods paper
for (( i = 0; i < ${#samples[@]}; i++ )); do
    sample=${samples[$i]}
    input_file="$data/trimmed/${sample}_R1_001.fastq.gz $data/trimmed/${sample}_R2_001.fastq.gz"
    output_dir=$results/$sample/alignments/star
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/$sample.bam
    if [[ ! -f $output_file ]]; then
        cpus=12
        tmp_file=$output_dir/${sample}_Aligned.out.sam
        runscript=star_${sample}.sh
        echo "STAR --runThreadN $cpus --genomeDir $stargenomedir --readFilesIn $input_file --readFilesCommand zcat --outFileNamePrefix $output_dir/${sample}_ --outFilterMultimapNmax 2 --sjdbGTFfile $hg19gtf" > $runscript
        echo "samtools view -ShuF4 $tmp_file | samtools sort -@ $cpus -m 16G - ${output_file/.bam}" >> $runscript
        echo "samtools index $output_file" >> $runscript
        echo "rm $tmp_file" >> $runscript
        bsub -J align -o align.%J.out -e align.%J.err -P $pid -R "select[mem>16] rusage[mem=16] span[hosts=1]" -n $cpus -K < $runscript &
    fi
done
wait


# remove reads spanning junction
# could use https://github.com/mjafin/splitNreads/blob/master/splitNreads.py
# rather than removing.
for (( i = 0; i < ${#samples[@]}; i++ )); do
    sample=${samples[$i]}
    input_file=$results/$sample/alignments/star/$sample.bam
    output_dir=$results/$sample/alignments/remove_n
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/$sample.bam
    if [[ ! -f $output_file ]]; then
        runscript=rm_n_${sample}.sh
        echo "samtools view -h $input_file | awk '\$6!~/N/' | samtools view -Sb - > $output_file" > $runscript
        bsub -J rm_n -o rm_n.%J.out -e rm_n.%J.err -P $pid -K < $runscript &
    fi
done
wait


# remove duplicate aligned reads
for (( i = 0; i < ${#samples[@]}; i++ )); do
    sample=${samples[$i]}
    input_file=$results/$sample/alignments/remove_n/$sample.bam
    output_dir=$results/$sample/alignments/remove_duplicates
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/$sample.bam
    if [[ ! -f $output_file ]]; then
        runscript=rmdup_${sample}.sh
        echo "samtools rmdup $input_file $output_file" > $runscript
        echo "samtools index $output_file" >> $runscript
        bsub -J rmdup -o rmdup.%J.out -e rmdup.%J.err -P $pid -K < $runscript &
    fi
done
wait


# call variants
for (( i = 0; i < ${#samples[@]}; i++ )); do
    sample=${samples[$i]}
    input_file=$results/$sample/alignments/remove_duplicates/$sample.bam
    output_dir=$results/$sample/variants/freebayes
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/$sample.vcf
    if [[ ! -f $output_file ]]; then
        cmd="freebayes -b $input_file -v $output_file -f $hg19reference -0"
        bsub -J fb -o fb.%J.out -e fb.%J.err -P $pid -K $cmd &
    fi
done
wait


# annotate with snpeff
for (( i = 0; i < ${#samples[@]}; i++ )); do
    sample=${samples[$i]}
    input_file=$results/$sample/variants/freebayes/$sample.vcf
    output_dir=$results/$sample/variants/snpeff
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/$sample.txt.gz
    if [[ ! -f $output_file ]]; then
        cmd="java -jar $snpeff eff -chr chr -minC 10 -no-downstream -no-intergenic -no-intron -no-upstream -noStats -v -c $snpeffconfig -o txt hg19 $input_file | gzip -c > $output_file"
        bsub -J snpeff -o snpeff.%J.out -e snpeff.%J.err -P $pid -K $cmd &
    fi
done
wait


# filter snpeff output
for (( i = 0; i < ${#samples[@]}; i++ )); do
    sample=${samples[$i]}
    input_file=$results/$sample/variants/snpeff/$sample.txt.gz
    output_dir=$results/$sample/variants/snpeff/filtered
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file=$output_dir/$sample.txt.gz
    if [[ ! -f $output_file ]]; then
        cmd="python ~/projects/ichiro_rnaseq/bin/gene_selection.py $input_file -e STOP_LOST_BUT_NOTHING_GAINED -e CODON_INSERTION -e EXON_DELETED -e FRAME_SHIFT -e NON_SYNONYMOUS_CODING -e NON_SYNONYMOUS_START -e NON_SYNONYMOUS_STOP -e START_LOST -e STOP_GAINED -e STOP_LOST -e NON_SYNONYMOUS_CODING | gzip -c > $output_file"
        bsub -J filtersnpeff -o fsnpeff.%J.out -e fsnpeff.%J.err -P $pid -K $cmd &
    fi
done
wait
