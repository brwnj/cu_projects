Peaktools
==============================================================================

prints None sometimes:
    
    chrY	None	None	MDX_22_AGTTCC_L003_R1_001_pos	None	+

Davidson; Human; T-Cell repertoire
==============================================================================

* SSAKE


Hits-Clip; hg18
==============================================================================

* ayb on index 51-56 only

* polya - cdc6; more to come from peter
* sum up the signal between annotated peaks
* rank between the comparison of the non-annotated regions - sum up the counts where not annotated, rank based on difference between case & control

* align; mask; align
* initial 22 bp of miRNA at 5' end
* align using novoalign
* of the alignments, lookup read in fastq, trim mapped 5' end that mapped
* remap new 3' piece
* associate read name with both sequences to give some idea of miRNA to mRNA binding

* do anything with trimming from the 3' end?

genetrix, genetricks - snRNA

* ayb with shorter and shorter (less priority)


Kamstock; mm9
==============================================================================

Next generation sequencing was conducted on triplicate RNA samples of 
nonmetastatic Dunn and the Fidler selected metastatic DLM8 subline. These 
samples were extracted from both cells grown in culture and tumors grown in 
mice. The samples were purified using dual rounds of OligodT column 
purification. In addition, DLM8 cells were treated in culture with 
2-deoxyglucose to determine the effect of this agent on cellular function.

Dunn13
Dunn1
Dunn2
Dunn3

DK15_Dunn_Tumor2
DK16_Dunn_Tumor3
DK17_Dunn_Tumor4

DK19_DLM8_Tumor1
DK20_DLM8_Tumor2
DK21_DLM8_Tumor3

DK23_Dunn_P9
DK25_Dunn_P10 *
DK13_Dunn_P18 *

DK1_DLM8
DK2_DLM8
DK27_DLM8

DK29_DLM8_2DG
DK31_DLM8_2DG
DK6_DLM8_2DG *

Novoalign
==============================================================================

making an index
----------------

download latest dbsnp for species
incorporate into reference

novoutil iupac dbsnp.vcf chr*.fasta.gz | gzip -c > hg19.dbsnp.fa.gz

create index

novoindex hg19.dbsnp.novoidx hg19.dbsnp.fa.gz


making an rnaseq index
-----------------------

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


Walter
==============================================================================

The big picture background is that we enrolled subjects in four groups:
Tuberculosis (TB)
Latent tuberculosis infection (LTBI)
Pneumonia (PNA)
No infection (No Infn)
 
We collected both whole blood (PAXgene) and PBMC. For preliminary data we
put 48 paired specimens from 48 subjects (12 from each group) on HuGene ST
1.1 arrays. These are the data we have available. This was paired in the sense 
that each individual subject had a PAXgene and PBMC specimen. The 48 PAXgene 
specimens are in an eset called PAXset.  The 48 PBMC specimens are in an eset 
called PBMCset.


















