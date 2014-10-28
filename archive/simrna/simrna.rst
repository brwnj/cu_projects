.. _simrna:

******************************************************************************
Simulated RNA
******************************************************************************

TODO
==============================================================================

* the remaining samples through rum

Tophat
==============================================================================

Script to queue::

    #! /usr/bin/env bash
    #BSUB -J tophat[1-8]
    #BSUB -e %J.%I.err
    #BSUB -o %J.%I.out
    #BSUB -n 12
    #BSUB -q normal
    #BSUB -R "span[hosts=1]"

    GENEMODEL=$HOME/projects/ref/hg19/Ensembl/GRCh37/Annotation/Archives/archive-2011-08-30-22-38-04/Genes/genes.gtf
    BOWTIEINDEX=$HOME/projects/ref/hg19/bowtie2/hg19

    # file format name is: human_1.fastq
    READS=$HOME/projects/simulatedrna/data/human_${LSB_JOBINDEX}.fastq

    # o = output dir
    # p = threads
    # g = only report single best alignment
    # coverage-search = most sensitive junction search
    # GTF = gene model reference
    # T = only align to gene model
    OPTIONS="-o ${LSB_JOBINDEX} -p 12 -g 1 --coverage-search --GTF $GENEMODEL -T"

    tophat2 $OPTIONS $BOWTIEINDEX $READS


Counts
==============================================================================

Script to queue::

    #! /usr/bin/env bash
    #BSUB -J genecounts[1-8]
    #BSUB -o %J.%I.out
    #BSUB -e %J.%I.err
    #BSUB -q normal
    #BSUB -w "exit(tophat[*], 0)"

    OPTIONS="-m intersection-nonempty -t exon -i gene_id"

    GENEMODEL=$HOME/projects/ref/hg19/Ensembl/GRCh37/Annotation/Genes/genes.chr.gtf

    BAM=$HOME/projects/simulatedrna/results/20120516/${LSB_JOBINDEX}/accepted_hits.bam
    OUTPUT=${LSB_JOBINDEX}.counts.txt

    samtools view $BAM | htseq-count $OPTIONS - $GENEMODEL > $OUTPUT


Running RNA-SeQC
==============================================================================

Script to queue::

    #! /usr/bin/env bash
    #BSUB -J rna_seqc
    #BSUB -e %J.err
    #BSUB -o %J.out
    #BSUB -q normal

    <<DOC
    Runs RNA-SeQC for the simulated reads.

    # verify contig names are consistent between bam, reference, and annotation
    # sort the bam by karotype using picard_tools/ReorderSam.jar
    # index bam and reference
    # create reference dictionary using picard_tools/CreateSequenceDictionary.jar
    DOC

    OUT=./simulated_rnaseq_qc/

    DATA=/vol1/home/brownj/projects/simulatedrna/data

    # .f contains 'chr'
    REFERENCE=$DATA/hg19.f.fa
    ANNOTATION=$DATA/hg19.f.gtf

    # make sample file
    SAMPLE_FILE=sample_desc.txt
    # header
    echo "SampleID"$'\t'"BamFile"$'\t'"Notes" > $SAMPLE_FILE
    for (( i = 1; i < 8; i++ )); do
        # tab delimited info
        echo -n "human_$i"$'\t' >> $SAMPLE_FILE
        echo -n "$DATA""/human_""$i"".karyosort.bam"$'\t' >> $SAMPLE_FILE
    	echo -n "simulated reads for human sample $i" >> $SAMPLE_FILE
    	# line break
    	echo "" >> $SAMPLE_FILE
    done

    OPTIONS="-n 5000 -o $OUT -r $REFERENCE -s $SAMPLE_FILE -singleEnd -t $ANNOTATION -ttype 2 -gld"

    java -jar ~/opt/bin/RNA-SeQC_v1.1.4.jar $OPTIONS

Test report: http://amc-einstein.ucdenver.pvt/~brownj/testReport/

Available metrics
------------------------------------------------------------------------------

* Read Metrics
    
    * Total, unique, duplicate reads
    * Alternative alignment reads
    * Read Length
    * Fragment Length mean and standard deviation
    * Read pairs: number aligned, unpaired reads, base mismatch rate for each pair mate, chimeric pairs
    * Vendor Failed Reads
    * Mapped reads and mapped unique reads
    * rRNA reads
    * Transcript-annotated reads (intragenic, intergenic, exonic, intronic)
    * Expression profiling efficiency (ratio of exon-derived reads to total reads sequenced)
    * Strand specificity

* Coverage

    * Mean coverage (reads per base)
    * Mean coefficient of variation
    * 5'/3' bias
    * Coverage gaps: count, length
    * Coverage Plots

* Downsampling
* GC Bias
* Correlation: 

    * Between sample(s) and a reference expression profile
    * When run with multiple samples, the correlation between every sample pair is reported


easyRNASeq
==============================================================================

::

    ## load the library
    library(easyRNASeq)

    ## load the chromosome sizes
    library(BSgenome.Hsapiens.UCSC.hg19)
    chr.sizes=as.list(seqlengths(Hsapiens))

    setwd("/Users/brownj/projects/simulated_rna")

    ## get the bam filenames
    bamfiles=dir(getwd(),pattern="*\\.bam$")

    ## run easyRNASeq to get an RNAseq object.
    rnaSeq <- easyRNASeq(filesDirectory=getwd(),
                         organism="Hsapiens",
                         chr.sizes=chr.sizes,
                         readLength=100L,
                         annotationMethod="biomaRt",
                         format="bam",
                         count="exons",
                         filenames=bamfiles[1],
                         outputFormat="RNAseq"
                         )

    gAnnot <- genomicAnnotation(rnaSeq)

    ## There are 244 "chromosomes" in that annotation, let's keep only what we need
    gAnnot <- gAnnot[space(gAnnot) %in% paste("chr",c(1:22,"X","Y","M"),sep=""),]
    save(gAnnot,file="gAnnot.rda")

    ## you could use it this way
    countTable <- easyRNASeq(filesDirectory=getwd(),
                             organism="Hsapiens",
                             chr.sizes=chr.sizes,
                             readLength=100L,
                             annotationMethod="rda",
                             annotationFile="gAnnot.rda",
                             format="bam",
                             count="genes",
                             summarization="geneModels",
                             filenames=bamfiles[1]
                             )

    ## applying DESeq
    ## Defining the conditions
    conditions <- c("control", "control","control","control","case","case","case","case")
    names(conditions) <- bamfiles

    ## running DESeq
    countDataSet <- easyRNASeq(filesDirectory=getwd(),
                               organism="Hsapiens",
                               chr.sizes=chr.sizes,
                               readLength=100L,
                               annotationMethod="env",
                               annotationObject=gAnnot,
                               format="bam",
                               count="genes",
                               summarization="geneModels",
                               filenames=bamfiles,
                               chr.sel="chr1",
                               outputFormat="DESeq",
                               normalize=TRUE,
                               conditions=conditions,
                               fitType="local",
                               method="per-condition"
                               )

    res = nbinomTest(countDataSet, "control", "case")
    write.table(res, file="simulated_rna.txt", sep="\t", quote=FALSE, row.names=FALSE)


Notes
==============================================================================

20120503
------------------------------------------------------------------------------

rna-seqc on bams
if bam with unique id, additional stats on read mapping locations


20120411
------------------------------------------------------------------------------

script needs to read fq, create list of uniq read ids?
from bam and gtf, use readid to see if read mapped to proper location
counts reads mapping outside of proper region
