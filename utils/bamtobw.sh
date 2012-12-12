#BSUB -J bamtobw
#BSUB -e bamtobw.%J.err
#BSUB -o bamtobw.%J.out

bedtools bamtobed -i test.bam | gzip -c > test.bed.gz
bedtools genomecov -strand + -bg -5 -ibam test.bam -g /path/to/sizes > test.pos.bedgraph
bedtools genomecov -strand - -bg -5 -ibam test.bam -g /path/to/sizes > test.neg.bedgraph
bedGraphToBigWig test.pos.bedgraph /path/to/sizes test.pos.bw
bedGraphToBigWig test.neg.bedgraph /path/to/sizes test.neg.bw
gzip -f test.*.bedgraph
