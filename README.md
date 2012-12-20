Peaktools
==============================================================================

prints None sometimes:
    
    chrY	None	None	MDX_22_AGTTCC_L003_R1_001_pos	None	+

Davidson; Human; T-Cell repertoire
==============================================================================

* SSAKE


Hits-Clip; hg18
==============================================================================

* polya - it looks like some sites are interesting and not showing in dexseq
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


Leinwand; mm9
==============================================================================

* per sample means appended to the DESeq output
* send more of the deseq output and qc including pca and clustering


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
DK25_Dunn_P10
DK13_Dunn_P18

DK1_DLM8
DK2_DLM8
DK27_DLM8

DK29_DLM8_2DG
DK31_DLM8_2DG
DK6_DLM8_2DG


chipseq chrhmm regions; run intervals based on regions specific to their region file;
plot p-vals for each group;