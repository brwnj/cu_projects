# source("http://bioconductor.org/biocLite.R")
# useDevel()
# biocLite("customProDB")
library(customProDB)

# this works fine but ensembl doesn't prepend "chr" to chromosome
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="feb2012.archive.ensembl.org", path="/biomart/martservice", archive=FALSE)
annotationpath <- "~/projects/davidson_transcriptome/data/ensembl"
PrepareAnnotationEnsembl(mart=ensembl, annotation_path=annotationpath, dbsnp="snp137")

# pepfasta <- "~/Downloads/protein.fasta"
# cdsfasta <- "~/Downloads/cds.fasta"
# annotationpath <- "~/projects/davidson_transcriptome/data/refseq"
# PrepareAnnotationRefseq(genome="hg19", cdsfasta, pepfasta, annotationpath, dbsnp="snp137")

path <- "/Users/brownj/projects/davidson_transcriptome/data/20131021/bams"
bams <- paste(path, "/", list.files(path, pattern="*bam$"), sep='')
load("~/projects/davidson_transcriptome/data/ensembl/ids.RData")
load("~/projects/davidson_transcriptome/data/ensembl/exon_anno.RData")
load("~/projects/davidson_transcriptome/data/ensembl/proseq.RData")
# load("~/projects/davidson_transcriptome/data/ensembl/dbsnpinCoding.RData")

rpkms <- sapply(bams, function(x) calculateRPKM(x, exon, proteincodingonly=TRUE, ids))
outfile <- "~/projects/davidson_transcriptome/data/20130927/rpkm_share.fasta"
pro <- OutputsharedPro(rpkms, cutoff=10, share_sample=2, proteinseq, outfile, ids)

path <- "/Users/brownj/projects/davidson_transcriptome/data/20131021/vcfs"
vcf_files <- paste(path, "/", list.files(path, pattern="*vcf$"), sep='')
vcfs <- lapply(vcf_files, function(x) InputVcf(x))
shared <- Multiple_VCF(vcfs, share_num=2)


bampath <- "~/projects/davidson_transcriptome/data/20131021/bams"
vcfpath <- "~/projects/davidson_transcriptome/data/20131021/vcfs"
annotationpath <- "~/projects/davidson_transcriptome/data/ensembl"
outfilepath <- "~/projects/davidson_transcriptome/data/20131021/"
outfilename <- "davidson"
easyRun_mul(bampath, vcfFile_path=vcfpath, annotation_path=annotationpath, 
            rpkm_cutoff=10, share_num=2, var_share_num=2, outfile_path=outfilepath,
            outfile_name=outfilename, INDEL=TRUE)