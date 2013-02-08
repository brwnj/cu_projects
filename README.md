#Davidson; Human; T-Cell repertoire

1alpha
2beta
3alpha
4beta
5alpha and 
6beta were spiked with some known sequence

#Hits-Clip; hg18

#Poly(A)

1. hybrids, run against mirbase. trim. mask. align to peaks. (what jay did)
2. after filtering to miRNA off, align to genome -r None, convert bam to fastq, realign to peak regions
3. run mirza on the output. figure out everything that jay just said...

#Walter

of the reads that mapped to tb, map those to human, call peaks, send bed for
each sample




#Novoalign

##making an index

download latest dbsnp for species
incorporate into reference

novoutil iupac dbsnp.vcf chr*.fasta.gz | gzip -c > hg19.dbsnp.fa.gz

create index

novoindex hg19.dbsnp.novoidx hg19.dbsnp.fa.gz

##making an rnaseq index

download refflat from ucsc

get rid of *random* and reorder columns

    bioawk -c header 'BEGIN{OFS=FS="\t"}\
        length($chrom)<6\
        {print $name2,$name,$chrom,$strand,$txStart,$txEnd,$cdsStart,$cdsEnd,$exonCount,$exonStarts,$exonEnds}' \
        refFlat_mm9_ensembl.txt.gz \
        | gzip -c > refFlat_mm9_ensembl.txt.fa.gz

make transcriptome using script in ~brownj/opt/bin

maketranscriptome fasta_dir refFlat read_length

create index for novoalign

novoindex mm9.mrna.novoidx refFlatRad45Num60kMin10Splices.fasta.gz refFlatRad45Num60kMin10Transcripts.fasta.gz

align

parse alignment
