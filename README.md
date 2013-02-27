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

#Tamim
* exome and genome to do cnv calls
* conifer (http://conifer.sourceforge.net/tutorial.html)
* 100-150 per year for 5 years; possibly whole genome
* early exome data didn't correlate with array data
* use existing exome data to see if we can find common cnvs using some software
* run conifer; see what it identifies

###some talk about pipelines
* james has been using galaxy and whatnot; he doesn't script so much
* we use unified genotyper, etc.

###array data
* process data from tamim's custom array
* database of cnv calls to identify common artifacts; reference: http://projects.tcag.ca/variation/

###moving forward
* 6 families, call variants, identify causal variant
* figure out what to do with conifer

#Walter

* of the reads that mapped to tb, map back to hg19 only the reads that fell in DE peaks
* peak counts for all samples
* reads counts for all samples
* average read length
* average width of peaks being tested in DESeq

* form consensus peak file from infected only; give number of total regions
* unafold to see nearby peaks form a stem-loop structure (http://www.idtdna.com/UNAFold)
* 54 bp; or 50-70 bp gap to account for stem
* locally align
* 3' utr, reverse strand alignment within transcriptome
* then take the gene hits into geo profiles to see if gene has been implicated in infectious disease
* lynne thinks 10 is a good number for candidates

blastn -query peaks.fa -strand minus -db dbpath -evalue 1000 -word_size 7 -gapopen 5 -gapextend 2 -reward 1 -penalty -3
blastn -query ~/projects/walter/data/20121124/peaks.fa -db ~/projects/walter/data/db/ensg_3utr -evalue 1000 -word_size 7 -gapopen 5 -gapextend 2 -reward 1 -penalty -3 -outfmt 6

chr1:1000273-1000313	hg18_ensGene_ENST00000321300	100.00	15	0	0	13	27	710	724	0.92	30.2
chr1:1000273-1000313	hg18_ensGene_ENST00000405804	100.00	15	0	0	22	36	1229	1243	0.92	30.2
chr1:1000273-1000313	hg18_ensGene_ENST00000400772	94.44	18	1	0	14	31	745	762	3.6	28.2
chr1:1000273-1000313	hg18_ensGene_ENST00000309965	100.00	14	0	0	16	29	129	116	3.6	28.2


#Novoalign

##making an index

download latest dbsnp for species
incorporate into reference

```
novoutil iupac dbsnp.vcf chr*.fasta.gz | gzip -c > hg19.dbsnp.fa.gz
```

create index

```
novoindex hg19.dbsnp.novoidx hg19.dbsnp.fa.gz
```

##making an rnaseq index

download refflat from ucsc

get rid of *random* and reorder columns

```
bioawk -c header 'BEGIN{OFS=FS="\t"}\
    length($chrom)<6\
    {print $name2,$name,$chrom,$strand,$txStart,$txEnd,$cdsStart,$cdsEnd,$exonCount,$exonStarts,$exonEnds}' \
    refFlat_mm9_ensembl.txt.gz \
    | gzip -c > refFlat_mm9_ensembl.txt.fa.gz
```

make transcriptome using script in ~brownj/opt/bin

```
maketranscriptome fasta_dir refFlat read_length
```

create index for novoalign

```
novoindex mm9.mrna.novoidx refFlatRad45Num60kMin10Splices.fasta.gz refFlatRad45Num60kMin10Transcripts.fasta.gz
```

align

parse alignment
