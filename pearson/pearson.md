#Pearson

need the fasta reference sorted by contig name
```
bioawk -c fastx '{print}' in.fa | sort -k1,1n | awk '{print ">"$1;print $2}'
```
then you need to make it a valid fasta
```
java -jar picard-tools-1.79/NormalizeFasta.jar I=tetrahymena_thermophila_sb210__mic__2_supercontigs.sorted.fasta  O=tetrahymena_thermophila_sb210__mic__2_supercontigs.sorted.normalized.fasta
```
index it
```
samtools faidx tetrahymena_thermophila_sb210__mic__2_supercontigs.sorted.normalized.fasta
```
create a dict of the fasta
```
java -jar picard-tools-1.79/CreateSequenceDictionary.jar R=tetrahymena_thermophila_sb210__mic__2_supercontigs.sorted.normalized.fasta GENOME_ASSEMBLY=tta_mic O=tetrahymena_thermophila_sb210__mic__2_supercontigs.sorted.normalized.fasta.dict
```
once complete
```
cd /vol1/home/brownj/projects/pearson/results/common
bedtools subtract -a DisA1_TGACCA_L001/DisA1_TGACCA_L001.mic.vcf -b B1868_ATCACG_L001/B1868_ATCACG_L001.mic.vcf > DisA1.uniq.mic.vcf
```
snpEff reference

+ add entry to snpEff.config
+ mv gtf to ~/ref/snpeff/<species>
+ mv the fasta annotation to genomes in ~/ref/snpeff/genomes
```
java -jar ~/opt/snpeff/snpEff.jar build -gtf22 -v -c ~/opt/snpeff/snpEff.config tta_mic
```
run snpEff
```
java -jar ~/opt/snpeff/snpEff.jar -i vcf -o txt -chr supercontig_ -minC 10 -minQ 30 -no-intergenic -ud 200 -hom -c ~/opt/snpeff/snpEff.config tta_mic DisA1.uniq.mic.vcf | awk '$13!="INTRON"' > DisA1_tta_mic.txt
```

##QC

+ http://amc-sandbox.ucdenver.edu/~brownj/tesla/_static/pearson/B1868_ATCACG_L001_R1_001_fastqc/fastqc_report.html
+ http://amc-sandbox.ucdenver.edu/~brownj/tesla/_static/pearson/B1868_ATCACG_L001_R2_001_fastqc/fastqc_report.html
+ http://amc-sandbox.ucdenver.edu/~brownj/tesla/_static/pearson/DisA1_TGACCA_L001_R1_001_fastqc/fastqc_report.html
+ http://amc-sandbox.ucdenver.edu/~brownj/tesla/_static/pearson/DisA1_TGACCA_L001_R2_001_fastqc/fastqc_report.html

##Alignment

+ Candidate list:

    + TTHERM_01308010 (Poc1)
    + TTHERM_01084190 (Bbs1)
    + TTHERM_00782070 (Bbs5)
    + TTHERM_00161270 (Bbc31)
    + TTHERM_00537420 (Fop1)
    + TTHERM_00216010 (Ftt18)
    + TTHERM_00068170 (Tcb2)

* reference genome (tetrahymena) and gene annotation exist
* discover unique mutations of the mutant where wildtype will represent the background
* b = wildtype, control, subtract what's different from the reference
* d = mutant, pool of f2s
* data resides in: pilot/120727
* fastqc, align, gatk
* sb210 strain
* list of genes of interest
* macro vs mic genome reference mapping
