#!/usr/bin/env bash
#BSUB -J bcell
#BSUB -e bcell.%J.err
#BSUB -o bcell.%J.out
#BSUB -q normal
#BSUB -n 1
#BSUB -P bennett

set -o nounset -o pipefail -o errexit -x
source /vol1/home/brownj/projects/bennett_nmo/bin/config.sh
trap cleanup EXIT


# always gzip everything upon completion
function cleanup () {
	echo "zipping fastqs"
	for f in `find $DATA -name *fastq`; do
		gzip -f $f
	done
}


SAMPLES=(control_ACATCG)

# filter quality
# $output_dir/${sample}_R[1-2]_quality-pass.fastq
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=filter_quality
    sample=${SAMPLES[$i]}

	input_file_1=$DATA/umi_added/${sample}_R1.fastq.gz
	input_file_2=$DATA/raw/${sample}_R2.fastq.gz

    output_dir=$DATA/quality_filtered
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi
    output_file_1=$output_dir/${sample}_R1_quality-pass.fastq.gz
	output_file_2=$output_dir/${sample}_R2_quality-pass.fastq.gz
	log_file_1=$RESULTS/logs/${sample}_R1_filter_quality.log
	log_file_2=$RESULTS/logs/${sample}_R2_filter_quality.log
    if [[ ! -f $output_file_1 ]]; then
        runscript=${jname}_${sample}_R1.sh
		echo "gunzip -f $input_file_1" > $runscript
		echo "FilterSeq.py quality -s ${input_file_1/.gz} -q $MINQUAL --nproc $NPROC --outdir $output_dir --outname ${sample}_R1 --log $log_file_1 --clean" >> $runscript
		echo "gzip -f ${input_file_1/.gz} ${output_file_1/.gz}" >> $runscript
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -n $NPROC -K < $runscript &
    fi
    if [[ ! -f $output_file_2 ]]; then
        runscript=${jname}_${sample}_R2.sh
		echo "gunzip -f $input_file_2" > $runscript
		echo "FilterSeq.py quality -s ${input_file_2/.gz} -q $MINQUAL --nproc $NPROC --outdir $output_dir --outname ${sample}_R2 --log $log_file_2 --clean" >> $runscript
		echo "gzip -f ${input_file_2/.gz} ${output_file_2/.gz}" >> $runscript
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -n $NPROC -K < $runscript &
    fi
done
wait


# identify primers and UMIs
# $output_dir/${sample}_R[1-2]_primers-pass.fastq
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
	jname=filter_primers
    sample=${SAMPLES[$i]}

	input_file_1=$DATA/quality_filtered/${sample}_R1_quality-pass.fastq.gz
	input_file_2=$DATA/quality_filtered/${sample}_R2_quality-pass.fastq.gz

    output_dir=$DATA/primer_filtered
    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi

    output_file_1=$output_dir/${sample}_R1_primers-pass.fastq.gz
	output_file_2=$output_dir/${sample}_R2_primers-pass.fastq.gz
	log_file_1=$RESULTS/logs/${sample}_R1_filter_primers.log
	log_file_2=$RESULTS/logs/${sample}_R2_filter_primers.log
    if [[ ! -f $output_file_1 ]]; then
        runscript=${jname}_${sample}_R1.sh
		echo "gunzip -f $input_file_1" > $runscript
		echo "MaskPrimers.py score -s ${input_file_1/.gz} -p $R1PRIMERS --mode cut --start $UMILENGTH --barcode --maxerror $R1_MAXERR --outdir $output_dir --outname ${sample}_R1 --nproc $NPROC --log $log_file_1 --clean" >> $runscript
		echo "gzip -f ${input_file_1/.gz} ${output_file_1/.gz}" >> $runscript
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -n $NPROC -K < $runscript &
    fi
    if [[ ! -f $output_file_2 ]]; then
        runscript=${jname}_${sample}_R2.sh
		echo "gunzip -f $input_file_2" > $runscript
		echo "MaskPrimers.py score -s ${input_file_2/.gz} -p $R2PRIMERS --mode cut --start $UMILENGTH --barcode --maxerror $R2_MAXERR --outdir $output_dir --outname ${sample}_R2 --nproc $NPROC --log $log_file_2 --clean" >> $runscript
		echo "gzip -f ${input_file_2/.gz} ${output_file_2/.gz}" >> $runscript
		bsub -J $jname -o $jname.%J.out -e $jname.%J.err -P $PI -n $NPROC -K < $runscript &
    fi
done
wait

exit


# Assign UIDs to read 1 sequences
echo "   3: PairSeq                $(date +'%H:%M %D')"
PairSeq.py -1 R1*primers-pass.fastq -2 R2*primers-pass.fastq -f BARCODE \
    --coord illumina --clean >> $RUNLOG

# Align sequence start positions by primer alignments
echo "   4: AlignSets              $(date +'%H:%M %D')"
#AlignSets.py table -p $R2_PRIMER_FILE --reverse --exec /usr/local/bin/muscle3.8.31_i86linux64 > /dev/null
AlignSets.py offset -s R2*pair-pass.fastq -d $R2_OFFSET_FILE --bf BARCODE --pf PRIMER \
    --nproc $NPROC >> $RUNLOG

# Build UID consensus sequences
echo "   5: BuildConsensus         $(date +'%H:%M %D')"
BuildConsensus.py -s R1*pair-pass.fastq --bf BARCODE --pf PRIMER --prcons $PRCONS \
    -q $MINQUAL --maxdiv $MAXDIV --nproc $NPROC --log ConsensusLogR1.log --clean >> $RUNLOG
BuildConsensus.py -s R2*pair-pass_align-pass.fastq --bf BARCODE --pf PRIMER \
    -q $MINQUAL --maxdiv $MAXDIV --nproc $NPROC --log ConsensusLogR2.log --clean >> $RUNLOG

# Assemble paired ends
echo "   6: AssemblePairs          $(date +'%H:%M %D')"
AssemblePairs.py align -1 R2*consensus-pass.fastq -2 R1*consensus-pass.fastq \
    --1f CONSCOUNT --2f PRCONS CONSCOUNT --coord presto --rc tail --maxerror $AP_MAXERR --alpha $ALPHA \
    --nproc $NPROC --log AssembleLog.log --outname Assembled --clean >> $RUNLOG

# Remove sequences with many Ns
echo "   7: FilterSeq missing      $(date +'%H:%M %D')"
FilterSeq.py missing -s Assembled_assemble-pass.fastq -n $FS_MISS --inner \
    --nproc $NPROC --log MissingLog.log >> $RUNLOG

# Rewrite header with minimum of CONSCOUNT
echo "   8: ParseHeaders collapse  $(date +'%H:%M %D')"
ParseHeaders.py collapse -s Assembled*missing-pass.fastq -f CONSCOUNT --act min --fasta > /dev/null

# Remove duplicate sequences
echo "   9: CollapseSeq            $(date +'%H:%M %D')"
CollapseSeq.py -s Assembled*reheader.fasta -n $CS_MISS --uf PRCONS \
    --cf CONSCOUNT --act sum --outname Assembled --inner >> $RUNLOG

# Filter to sequences with at least 2 supporting sources
echo "  10: SplitSeq group         $(date +'%H:%M %D')"
SplitSeq.py group -s Assembled_collapse-unique.fasta -f CONSCOUNT --num 2 >> $RUNLOG

# Process log files
echo "  11: ParseLog               $(date +'%H:%M %D')"
ParseLog.py -l QualityLogR[1-2].log -f ID QUALITY > /dev/null &
ParseLog.py -l PrimerLogR[1-2].log -f ID BARCODE PRIMER ERROR > /dev/null &
ParseLog.py -l ConsensusLogR[1-2].log -f BARCODE BCCOUNT CONSCOUNT PRIMER PRCOUNT PRFREQ DIVERSITY > /dev/null &
ParseLog.py -l AssembleLog.log -f ID OVERLAP LENGTH PVAL ERROR > /dev/null &
ParseLog.py -l MissingLog.log -f ID MISSING > /dev/null &
ParseLog.py -l Pipeline.log -f END SEQUENCES PAIRS SETS PASS FAIL UNIQUE DUPLICATE UNDETERMINED PARTS OUTPUT > /dev/null &
wait

# Zip intermediate and log files
if $ZIP_FILES; then
    tar -cf LogFiles.tar *.log
    gzip LogFiles.tar
    rm *.log

    tar -cf TempFiles.tar R[1-2]_*.fastq *under* *duplicate* *undetermined* *reheader*
    gzip TempFiles.tar
    rm R[1-2]_*.fastq *under* *duplicate* *undetermined* *reheader*
fi

# End
echo -e "DONE\n"
cd ../

