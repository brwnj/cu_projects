setwd("~/projects/polya/data/02052013")
for(i in 1:6){
    t <- read.table(file=paste0(i,".tab"), header=FALSE, sep="\t")
    f <- read.table(file=paste0(i,".f.tab"), header=FALSE, sep="\t")
    pdf(file=paste0(i,".pdf"))
    plot(t$V1, t$V2, type="l", xlab="Distance (mm)", ylab="Absorbance", col="red")
    abline(v=f$V1, col="lightgray", lty=3)
    text(f$V1, max(t$V2), labels=seq_along(f$V1), col="lightgray")
    dev.off()
}