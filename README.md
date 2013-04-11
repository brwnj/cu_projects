#Davidson; Human; T-Cell repertoire

+ awaiting reply email

#Hits-Clip; hg18

PK61-PK63 the last run
++ still need to be run!

#Poly(A); hg18

+ alias bkill
+ classifying peaks
+ run fisher test on new sample pairs

### finding novel poly(a) sites

+ est libraries that polya primed
+ high confidence peaks at the 3' end of est as another class in addition to the other way
+ look near the annotated peak for new polya site

* see paper from jay
* compare for each site then for each gene
* a[a,t]taaa
if 8 or more As follow the 3 As, class 2
class one has fewer than 8 As after that stretch
* classifying the other groups (3, 4) -- call peaks and classify

category
1   has A[A,T]TAAA; NOT A-rich downstream from this cleavage site
2   has A[A,T]TAAA; with A-rich sequence downstream
3   lacks A[A,T]TAAA; no A-rich region downstream
4   only downstream A-rich sequence

100 bp from known polya sites

filter out peaks with less than 10 reads of support
canonical PAS should be found in -10 to -30 of the cleavage site
no more than three canonical PAS should be located in the upstream window

#Tamim
* conifer on different samples and compare to microarray data

#Walter

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

#A5

had to remove check for phred64, tis not working
replaced idba with idba_ud in ngopt/bin
