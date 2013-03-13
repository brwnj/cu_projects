#Davidson; Human; T-Cell repertoire

1alpha
2beta
3alpha
4beta
5alpha and 
6beta were spiked with some known sequence

sameple:total reads
1:51307090
2:45512106
3:38216696
4:61373826
5:48097530
6:42528694

#Hits-Clip; hg18
* peter -- correlate microarray abundance to miRNA abundance levels
    predicted network
    miRNA9 high miRNA193 low
    targets appear to be regulated by a given miRNA
    do we see genes going up, down, etc...
    erp+
    contigency table
            9   193
    dead
    alive
    chi-square test
    contigency table
                    miRNA_high  miRNA_low
    dead_outcome
    alive_outcome
    
    * define which miRNAs are implicated in regulating estrogen and tamoxiphen
    PK63 is tamoxiphen resistant
* continue getting PK samples through the pipeline

#Poly(A); hg18
* with and without umi mapped reads counts
* for each gene, contigency table to run fishers exact test
            gene p1 p2 p3 p4
c1
c2
c3

* testing gene level changes -- one condition versus one other
* testing counts with and without umi -- does it make a difference?
* can we run something faster like bowtie and get the same result as novoalign

### finding novel poly(a) sites
* see paper from jay
* compare for each site then for each gene
* a[a,t]taaa
if 8 or more As follow the 3 As, class 2
class one has fewer than 8 As after that stretch
* classifying the other groups (3, 4) -- call peaks and classify

could just call peaks, remove those within 50 bp from known polya sites
filter out peaks without characterization or around long stretches of As
doesn't get you all possible peaks though...

#Tamim
* conifer on different samples and compare to microarray data

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
