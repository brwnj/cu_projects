#!/usr/bin/env bash
#BSUB -J bcell
#BSUB -e bcell.%J.err
#BSUB -o bcell.%J.out
#BSUB -q normal
#BSUB -n 1
#BSUB -P bennett

set -o nounset -o pipefail -o errexit -x
source /vol1/home/brownj/projects/bennett_nmo/bin/config.sh

<<DOC
+ when sequences come from the core, transfer and rename to sample_R?.fastq.gz
+ prepend the UMI from the index read onto R1 using prepend_umi.py; reads go into umi_added folder
+ add sample into config.sh SAMPLES array
+ run presto.sh

there's a lot of overhead due to gunzip/gzip actions. ideally,
you'd combine all of the actions into a single job, gzipping
once, but i was testing each step individually.

DOC

# always gzip everything upon completion
function cleanup () {
	echo "zipping fastqs"
	for f in `find $DATA/*/* -name *fast[qa]`; do
		echo "gzip -f $f" | bsez -J cleaningup -P $PI
	done
}
trap cleanup EXIT


# filter quality
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=quality_pass
    sample=${SAMPLES[$i]}

	pipeline_log=$RESULTS/logs/${sample}_pipeline.log

	input_file_1=$DATA/umi_added/${sample}_R1.fastq.gz
	input_file_2=$DATA/raw/${sample}_R2.fastq.gz

	if [[ -f $input_file_1 && -f $input_file_2 ]]; then

	    output_dir=$DATA/$jname
	    if [[ ! -d $output_dir ]]; then
	        mkdir -p $output_dir
	    fi
	    output_file_1=$output_dir/${sample}_R1_quality-pass.fastq.gz
		output_file_2=$output_dir/${sample}_R2_quality-pass.fastq.gz
		log_file_1=$RESULTS/logs/${sample}_R1_${jname}.log
		log_file_2=$RESULTS/logs/${sample}_R2_${jname}.log
	    if [[ ! -f $output_file_1 ]]; then
	        runscript=${jname}_${sample}_R1.sh
			echo "gunzip -f $input_file_1" > $runscript
			echo "FilterSeq.py quality -s ${input_file_1/.gz} -q $MINQUAL --nproc $NPROC --outdir $output_dir --outname ${sample}_R1 --log $log_file_1 --clean >> ${pipeline_log}" >> $runscript
			echo "gzip -f ${input_file_1/.gz} ${output_file_1/.gz}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -n $NPROC -K -R "span[hosts=1]" < $runscript &
	    fi
	    if [[ ! -f $output_file_2 ]]; then
	        runscript=${jname}_${sample}_R2.sh
			echo "gunzip -f $input_file_2" > $runscript
			echo "FilterSeq.py quality -s ${input_file_2/.gz} -q $MINQUAL --nproc $NPROC --outdir $output_dir --outname ${sample}_R2 --log $log_file_2 --clean >> ${pipeline_log}" >> $runscript
			echo "gzip -f ${input_file_2/.gz} ${output_file_2/.gz}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -n $NPROC -K -R "span[hosts=1]" < $runscript &
	    fi

	fi
done
wait


# identify primers and UMIs
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=primers_pass
    sample=${SAMPLES[$i]}

	pipeline_log=$RESULTS/logs/${sample}_pipeline.log

	input_file_1=$DATA/quality_pass/${sample}_R1_quality-pass.fastq.gz
	input_file_2=$DATA/quality_pass/${sample}_R2_quality-pass.fastq.gz

	if [[ -f $input_file_1 && $input_file_2 ]]; then

	    output_dir=$DATA/$jname
	    if [[ ! -d $output_dir ]]; then
	        mkdir -p $output_dir
	    fi

	    output_file_1=$output_dir/${sample}_R1_primers-pass.fastq.gz
		output_file_2=$output_dir/${sample}_R2_primers-pass.fastq.gz
		log_file_1=$RESULTS/logs/${sample}_R1_${jname}.log
		log_file_2=$RESULTS/logs/${sample}_R2_${jname}.log
	    if [[ ! -f $output_file_1 ]]; then
	        runscript=${jname}_${sample}_R1.sh
			echo "gunzip -f $input_file_1" > $runscript
			echo "MaskPrimers.py score -s ${input_file_1/.gz} -p $R1PRIMERS --mode cut --start $UMILENGTH --barcode --maxerror $R1_MAXERR --outdir $output_dir --outname ${sample}_R1 --nproc $NPROC --log $log_file_1 --clean >> ${pipeline_log}" >> $runscript
			echo "gzip -f ${input_file_1/.gz} ${output_file_1/.gz}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -n $NPROC -K -R "span[hosts=1]" < $runscript &
	    fi
	    if [[ ! -f $output_file_2 ]]; then
	        runscript=${jname}_${sample}_R2.sh
			echo "gunzip -f $input_file_2" > $runscript
			echo "MaskPrimers.py score -s ${input_file_2/.gz} -p $R2PRIMERS --mode cut --start 0 --maxerror $R2_MAXERR --outdir $output_dir --outname ${sample}_R2 --nproc $NPROC --log $log_file_2 --clean >> ${pipeline_log}" >> $runscript
			echo "gzip -f ${input_file_2/.gz} ${output_file_2/.gz}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -n $NPROC -K -R "span[hosts=1]" < $runscript &
	    fi

	fi
done
wait


# pair the sequences based on UMI
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=pair_pass
    sample=${SAMPLES[$i]}

	pipeline_log=$RESULTS/logs/${sample}_pipeline.log

	input_file_1=$DATA/primers_pass/${sample}_R1_primers-pass.fastq.gz
	input_file_2=$DATA/primers_pass/${sample}_R2_primers-pass.fastq.gz

	if [[ -f $input_file_1 && $input_file_2 ]]; then

	    output_dir=$DATA/$jname
	    if [[ ! -d $output_dir ]]; then
	        mkdir -p $output_dir
	    fi

	    output_file_1=$output_dir/${sample}_R1_pair-pass.fastq.gz
		output_file_2=$output_dir/${sample}_R2_pair-pass.fastq.gz
		log_file_1=$RESULTS/logs/${sample}_R1_${jname}.log
		log_file_2=$RESULTS/logs/${sample}_R2_${jname}.log
	    if [[ ! -f $output_file_1 ]]; then
	        runscript=${jname}_${sample}.sh
			echo "set -o nounset -o pipefail -o errexit -x" > $runscript
			echo "gunzip -f $input_file_1 $input_file_2" >> $runscript
			echo "PairSeq.py -1 ${input_file_1/.gz} -2 ${input_file_2/.gz} -f BARCODE --coord illumina --clean --outdir $output_dir >> ${pipeline_log}" >> $runscript
			# can't specify outname due to multiple files being passed
			echo "mv ${output_file_1/_pair-pass.fastq.gz/_primers-pass_pair-pass.fastq} ${output_file_1/.gz}" >> $runscript
			echo "mv ${output_file_2/_pair-pass.fastq.gz/_primers-pass_pair-pass.fastq} ${output_file_2/.gz}" >> $runscript
			echo "gzip -f ${input_file_1/.gz} ${input_file_2/.gz} ${output_file_1/.gz} ${output_file_2/.gz}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K < $runscript &
	    fi

	fi
done
wait


# create offset table for R2 primers
if [[ ! -f $R2PRIMEROFFSET ]]; then
	AlignSets.py table --outdir $DATA --outname ${R2PRIMEROFFSET/_offsets-reverse.tab} -p $R2PRIMERS --reverse
fi


# align sequence start positions by primer alignments
# muscle path is hardcoded in AlignSets.py...
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=align_pass
    sample=${SAMPLES[$i]}

	pipeline_log=$RESULTS/logs/${sample}_pipeline.log

	input_file=$DATA/pair_pass/${sample}_R2_pair-pass.fastq.gz

	if [[ -f $input_file ]]; then

	    output_dir=$DATA/$jname
	    if [[ ! -d $output_dir ]]; then
	        mkdir -p $output_dir
	    fi

	    output_file=$output_dir/${sample}_R2_align-pass.fastq.gz
		log_file=$RESULTS/logs/${sample}_R2_${jname}.log
	    if [[ ! -f $output_file ]]; then
	        runscript=${jname}_${sample}.sh
			echo "set -o nounset -o pipefail -o errexit -x" > $runscript
			echo "gunzip -f $input_file" >> $runscript
			echo "AlignSets.py offset -s ${input_file/.gz} -d $R2PRIMEROFFSET --bf BARCODE --pf PRIMER --outdir $output_dir --outname ${sample}_R2 --nproc $NPROC --log $log_file --clean >> ${pipeline_log}" >> $runscript
			echo "gzip -f ${input_file/.gz} ${output_file/.gz}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -n $NPROC -K -R "span[hosts=1]" < $runscript &
	    fi

	fi
done
wait


# build UMI consensus
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=consensus_pass
    sample=${SAMPLES[$i]}

	pipeline_log=$RESULTS/logs/${sample}_pipeline.log

	input_file_1=$DATA/pair_pass/${sample}_R1_pair-pass.fastq.gz
	input_file_2=$DATA/align_pass/${sample}_R2_align-pass.fastq.gz

	if [[ -f $input_file_1 && $input_file_2 ]]; then

	    output_dir=$DATA/$jname
	    if [[ ! -d $output_dir ]]; then
	        mkdir -p $output_dir
	    fi

	    output_file_1=$output_dir/${sample}_R1_consensus-pass.fastq.gz
		output_file_2=$output_dir/${sample}_R2_consensus-pass.fastq.gz
		log_file_1=$RESULTS/logs/${sample}_R1_${jname}.log
		log_file_2=$RESULTS/logs/${sample}_R2_${jname}.log
	    if [[ ! -f $output_file_1 ]]; then
	        runscript=${jname}_${sample}_R1.sh
			echo "gunzip -f $input_file_1" > $runscript
			echo "BuildConsensus.py -s ${input_file_1/.gz} --bf BARCODE --pf PRIMER --prcons $PRCONS -q $MINQUAL --maxdiv $MAXDIV --nproc $NPROC --outdir $output_dir --outname ${sample}_R1 --log $log_file_1 --clean >> ${pipeline_log}" >> $runscript
			echo "gzip -f ${input_file_1/.gz} ${output_file_1/.gz}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -n $NPROC -K -R "span[hosts=1]" < $runscript &
	    fi
	    if [[ ! -f $output_file_2 ]]; then
	        runscript=${jname}_${sample}_R2.sh
			echo "gunzip -f $input_file_2" > $runscript
			echo "BuildConsensus.py -s ${input_file_2/.gz} --bf BARCODE --pf PRIMER -q $MINQUAL --maxdiv $MAXDIV --nproc $NPROC --outdir $output_dir --outname ${sample}_R2 --log $log_file_2 --clean >> ${pipeline_log}" >> $runscript
			echo "gzip -f ${input_file_2/.gz} ${output_file_2/.gz}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -n $NPROC -K -R "span[hosts=1]" < $runscript &
	    fi

	fi
done
wait


# Assemble paired ends
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=assemble_pass
    sample=${SAMPLES[$i]}

	pipeline_log=$RESULTS/logs/${sample}_pipeline.log

	input_file_1=$DATA/consensus_pass/${sample}_R1_consensus-pass.fastq.gz
	input_file_2=$DATA/consensus_pass/${sample}_R2_consensus-pass.fastq.gz

	if [[ -f $input_file_1 && $input_file_2 ]]; then

	    output_dir=$DATA/$jname
	    if [[ ! -d $output_dir ]]; then
	        mkdir -p $output_dir
	    fi

	    output_file=$output_dir/${sample}_assemble-pass.fastq.gz
		log_file=$RESULTS/logs/${sample}_${jname}.log
	    if [[ ! -f $output_file ]]; then
	        runscript=${jname}_${sample}.sh
			echo "gunzip -f $input_file_1 $input_file_2" > $runscript
			echo "AssemblePairs.py align -1 ${input_file_2/.gz} -2 ${input_file_1/.gz} --1f CONSCOUNT --2f PRCONS CONSCOUNT --coord presto --rc tail --maxerror $AP_MAXERR --alpha $ALPHA --nproc $NPROC --outdir $output_dir --outname ${sample} --log $log_file --clean >> ${pipeline_log}" >> $runscript
			echo "gzip -f ${input_file_1/.gz} ${input_file_2/.gz} ${output_file/.gz}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -n $NPROC -K -R "span[hosts=1]" < $runscript &
	    fi

	fi
done
wait


# Remove sequences with many Ns
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=missing_pass
    sample=${SAMPLES[$i]}

	pipeline_log=$RESULTS/logs/${sample}_pipeline.log

	input_file=$DATA/assemble_pass/${sample}_assemble-pass.fastq.gz

	if [[ -f $input_file ]]; then

	    output_dir=$DATA/$jname
	    if [[ ! -d $output_dir ]]; then
	        mkdir -p $output_dir
	    fi

	    output_file=$output_dir/${sample}_missing-pass.fastq.gz
		log_file=$RESULTS/logs/${sample}_${jname}.log
	    if [[ ! -f $output_file ]]; then
	        runscript=${jname}_${sample}.sh
			echo "gunzip -f $input_file" > $runscript
			echo "FilterSeq.py missing -s ${input_file/.gz} -n $FS_MISS --inner --nproc $NPROC --outdir $output_dir --outname ${sample} --log $log_file --clean >> ${pipeline_log}" >> $runscript
			echo "gzip -f ${input_file/.gz}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -n $NPROC -K -R "span[hosts=1]" < $runscript &
	    fi

	fi
done
wait


# Rewrite header with minimum of CONSCOUNT
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=reheader
    sample=${SAMPLES[$i]}

	pipeline_log=$RESULTS/logs/${sample}_pipeline.log

	input_file=$DATA/missing_pass/${sample}_missing-pass.fastq.gz

	if [[ -f $input_file ]]; then

	    output_dir=$DATA/$jname
	    if [[ ! -d $output_dir ]]; then
	        mkdir -p $output_dir
	    fi

	    output_file=$output_dir/${sample}_reheader.fasta.gz
		log_file=$RESULTS/logs/${sample}_${jname}.log
	    if [[ ! -f $output_file ]]; then
	        runscript=${jname}_${sample}.sh
			echo "gunzip -f $input_file" > $runscript
			echo "ParseHeaders.py collapse -s ${input_file/.gz} -f CONSCOUNT --act min --fasta --outdir $output_dir --outname $sample --clean >> ${pipeline_log}" >> $runscript
			echo "gzip -f ${input_file/.gz} ${output_file/.gz}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K < $runscript &
	    fi

	fi
done
wait


# Remove duplicate sequences
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=collapse_unique
    sample=${SAMPLES[$i]}

	pipeline_log=$RESULTS/logs/${sample}_pipeline.log

	input_file=$DATA/reheader/${sample}_reheader.fasta.gz

	if [[ -f $input_file ]]; then

	    output_dir=$DATA/$jname
	    if [[ ! -d $output_dir ]]; then
	        mkdir -p $output_dir
	    fi

	    output_file=$output_dir/${sample}_collapse-unique.fasta.gz
		log_file=$RESULTS/logs/${sample}_${jname}.log
	    if [[ ! -f $output_file ]]; then
	        runscript=${jname}_${sample}.sh
			echo "gunzip -f $input_file" > $runscript
			echo "CollapseSeq.py -s ${input_file/.gz} -n $CS_MISS --uf PRCONS --cf CONSCOUNT --act sum --inner --outdir $output_dir --outname $sample --log $log_file --clean >> ${pipeline_log}" >> $runscript
			echo "gzip -f ${input_file/.gz} ${output_file/.gz}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K < $runscript &
	    fi

	fi
done
wait


# Keep sequences with at least 2 supporting sources
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=at_least_2
    sample=${SAMPLES[$i]}

	pipeline_log=$RESULTS/logs/${sample}_pipeline.log

	input_file=$DATA/collapse_unique/${sample}_collapse-unique.fasta.gz

	if [[ -f $input_file ]]; then

	    output_dir=$DATA/$jname
	    if [[ ! -d $output_dir ]]; then
	        mkdir -p $output_dir
	    fi

	    output_file=$output_dir/${sample}_atleast-2.fasta.gz
		log_file=$RESULTS/logs/${sample}_${jname}.log
	    if [[ ! -f $output_file ]]; then
	        runscript=${jname}_${sample}.sh
			echo "gunzip -f $input_file" > $runscript
			echo "SplitSeq.py group -s ${input_file/.gz} -f CONSCOUNT --num 2 --outdir $output_dir --outname $sample --clean >> ${pipeline_log}" >> $runscript
			echo "gzip -f ${input_file/.gz} ${output_file/.gz}" >> $runscript
			bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K < $runscript &
	    fi

	fi
done
wait


# Process log files
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=tables
    sample=${SAMPLES[$i]}
	outdir=$RESULTS/$sample

	if [[ ! -d $outdir ]]; then
		mkdir -p $outdir
	fi

	input_file=$RESULTS/logs/${sample}_R1_quality_pass.log
	output_file=$outdir/$sample/${sample}_R1_quality_pass_table.tab
	if [[ -f $input_file && ! -f $output_file ]]; then
		cmd="ParseLog.py -l $input_file -f ID QUALITY --outdir $outdir"
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
	fi

	input_file=$RESULTS/logs/${sample}_R2_quality_pass.log
	output_file=$outdir/$sample/${sample}_R2_quality_pass_table.tab
	if [[ -f $input_file && ! -f $output_file ]]; then
		cmd="ParseLog.py -l $input_file -f ID QUALITY --outdir $outdir"
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
	fi

	input_file=$RESULTS/logs/${sample}_R1_primers_pass.log
	output_file=$outdir/$sample/${sample}_R1_primers_pass_table.tab
	if [[ -f $input_file && ! -f $output_file ]]; then
		cmd="ParseLog.py -l $input_file -f ID BARCODE PRIMER ERROR --outdir $outdir"
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
	fi

	input_file=$RESULTS/logs/${sample}_R2_primers_pass.log
	output_file=$outdir/$sample/${sample}_R2_primers_pass_table.tab
	if [[ -f $input_file && ! -f $output_file ]]; then
		cmd="ParseLog.py -l $input_file -f ID BARCODE PRIMER ERROR --outdir $outdir"
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
	fi

	input_file=$RESULTS/logs/${sample}_R1_consensus_pass.log
	output_file=$outdir/$sample/${sample}_R1_consensus_pass_table.tab
	if [[ -f $input_file && ! -f $output_file ]]; then
		cmd="ParseLog.py -l $input_file -f BARCODE BCCOUNT CONSCOUNT PRIMER PRCOUNT PRFREQ DIVERSITY --outdir $outdir"
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
	fi

	input_file=$RESULTS/logs/${sample}_R2_consensus_pass.log
	output_file=$outdir/$sample/${sample}_R2_consensus_pass_table.tab
	if [[ -f $input_file && ! -f $output_file ]]; then
		cmd="ParseLog.py -l $input_file -f BARCODE BCCOUNT CONSCOUNT PRIMER PRCOUNT PRFREQ DIVERSITY --outdir $outdir"
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
	fi

	input_file=$RESULTS/logs/${sample}_assemble_pass.log
	output_file=$outdir/$sample/${sample}_assemble_pass_table.tab
	if [[ -f $input_file && ! -f $output_file ]]; then
		cmd="ParseLog.py -l $input_file -f ID OVERLAP LENGTH PVAL ERROR --outdir $outdir"
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
	fi

	input_file=$RESULTS/logs/${sample}_missing_pass.log
	output_file=$outdir/$sample/${sample}_missing_pass_table.tab
	if [[ -f $input_file && ! -f $output_file ]]; then
		cmd="ParseLog.py -l $input_file -f ID MISSING --outdir $outdir"
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
	fi

	input_file=$RESULTS/logs/${sample}_pipeline.log
	output_file=$outdir/$sample/${sample}_pipeline_table.tab
	if [[ -f $input_file && ! -f $output_file ]]; then
		cmd="ParseLog.py -l $input_file -f END SEQUENCES PAIRS SETS PASS FAIL UNIQUE DUPLICATE UNDETERMINED PARTS OUTPUT --outdir $outdir"
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K $cmd &
	fi

	input_file=$DATA/at_least_2/${sample}_atleast-2.fasta.gz
	output_file=$outdir/$sample/${sample}_atleast-2_header.tab
	if [[ -f $input_file && ! -f $output_file ]]; then
        runscript=${jname}_${sample}.sh
		echo "gunzip -f $input_file" > $runscript
		echo "ParseHeaders.py table -s ${input_file/.gz} -f ID PRCONS CONSCOUNT DUPCOUNT --outdir $outdir" >> $runscript
		echo "gzip -f ${input_file/.gz}" >> $runscript
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -K < $runscript &
	fi

done
wait

# gzip all of the logs
# for f in $RESULTS/logs/*.log; do
# 	echo "gzip -f $f" | bsez -J cleaningup -P $PI
# done
#
# # gzip the summary tables
# for f in $RESULTS/*/*.tab; do
# 	echo "gzip -f $f" | bsez -J cleaningup -P $PI
# done
