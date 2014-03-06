
usage <- "Usage: Rscript bam targetsbed name"

args <- commandArgs(TRUE)
if (length(args) == 0)
    stop(usage)

library(TEQC)

readsfile = args[1]
targetsfile = args[2]
sample = args[3]

reads = get.reads(readsfile, filetype="bam")
targets = get.targets(targetsfile, chrcol=1, startcol=2, endcol=3, skip=0)

Coverage = coverage.target(reads, targets, perTarget=TRUE, perBase=TRUE)

pdf(paste0(sample, "_coverage.pdf"))
coverage.hist(Coverage$coverageTarget, covthreshold=10)
dev.off()
