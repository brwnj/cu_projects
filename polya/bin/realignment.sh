#! /usr/bin/env bash
#BSUB -J bowtie
#BSUB -e bowtie.%J.err
#BSUB -o bowtie.%J.out
#BSUB -q normal
#BSUB -R "select[mem>12] rusage[mem=12] span[hosts=1]"
#BSUB -n 1

<<DOC
Align using Bowtie, suppressing all reads that align more than once.
DOC

set -o nounset -o pipefail -o errexit -x

bidx=/vol3/home/jhessel/projects/bowtie/indices/hg18
mpg1=/vol1/home/brownj/projects/polya/data/20121008/MPG1.ACTG.fastq
mp41=/vol1/home/brownj/projects/polya/data/20121206/MP41.fastq

results=/vol1/home/brownj/projects/polya/results/common
mpg1bam=$results/MPG1.ACTG/MPG1.uniques.bam
mp41bam=$results/MP41/MP41.uniques.bam

# gunzip $mp41
# gunzip $mpg1bam
# bowtie -p4 -q -m1 --sam $bidx $mpg1 | samtools view -ShuF4 - | samtools sort -o - temp -m 9500000000 > $mpg1bam
bowtie -p4 -q -m1 --sam $bidx $mp41 | samtools view -ShuF4 - | samtools sort -o - temp -m 9500000000 > $mp41bam
gzip $mpg1 $mp41

sizes=/vol1/home/brownj/ref/hg18/hg18.sizes

bam2bw $mpg1bam $sizes
bam2bw $mp41bam $sizes

# track type=bigWig name='MP41 neg' description='MP41 neg' bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/polya/20121207/MP41.uniques.neg.bw maxHeightPixels=15:50:35 color=228,26,28 visibility=full
# track type=bigWig name='MP41 pos' description='MP41 pos' bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/polya/20121207/MP41.uniques.pos.bw maxHeightPixels=15:50:35 color=228,26,28 visibility=full
# track type=bigWig name='MPG1.ACTG neg' description='MPG1 neg' bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/polya/20121207/MPG1.uniques.neg.bw maxHeightPixels=15:50:35 color=55,126,184 visibility=full
# track type=bigWig name='MPG1.ACTG pos' description='MPG1 pos' bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/polya/20121207/MPG1.uniques.pos.bw maxHeightPixels=15:50:35 color=55,126,184 visibility=full