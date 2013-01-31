Davidson; Human; T-Cell repertoire
==============================================================================

* SSAKE


Hits-Clip; hg18
==============================================================================

* polya - cdc6; more to come from peter
* sum up the signal between annotated peaks
* rank between the comparison of the non-annotated regions - sum up the counts where not annotated, rank based on difference between case & control

* align; mask; align
* initial 22 bp of miRNA at 5' end
* align using novoalign
* of the alignments, lookup read in fastq, trim mapped 5' end that mapped
* remap new 3' piece
* associate read name with both sequences to give some idea of miRNA to mRNA binding


Kamstock; rnaseq; mm9
==============================================================================

* analyze cuffdiff results


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

align

parse alignment


Walter
==============================================================================

of the reads that mapped to tb, map those to human, call peaks, send bed for
each sample