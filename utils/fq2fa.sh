#!/usr/bin/env bash
set -o nounset -o errexit
function usage () {
    echo "convert fastq(.gz) to gzipped fasta. submits to normal queue.
    "
    echo "Usage: $(basename $0) <fastq> <project_id>"
    exit 1
}

if [ $# -lt 2 ]; then
    usage
fi

fastq=$1
projectid=$2
out=${fastq%.f*}.fa.gz
runscript="fastq2fasta.sh"
cat <<runscript >$runscript
#BSUB -J fq2fa
#BSUB -e fq2fa.%J.err
#BSUB -o fq2fa.%J.out
#BSUB -P $projectid

bioawk -c fastx '{print ">"\$name"\n"\$seq}' $fastq | gzip -c > $out
runscript

job=$(bsub < $runscript)
jobid=$(echo $job | sed -rn 's/.*<([0-9]+)>.*/\1/p')

mv fastq2fasta.sh fastq2fasta.$jobid.sh
echo $jobid