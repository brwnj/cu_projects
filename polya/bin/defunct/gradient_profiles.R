setwd("~/projects/polya/data/02052013/gradients/")
samples=c("ctrl-2-1", "ctrl-2-2", "ctrl-2-3", "msi-9-1", "msi-9-2", "msi-9-3")
for(i in seq_along(samples)){
    t <- read.table(file=paste0(samples[i],".tab"), header=FALSE, sep="\t")
    f <- read.table(file=paste0(samples[i],".f.tab"), header=FALSE, sep="\t")
    pdf(file=paste0(samples[i],".pdf"))
    plot(t$V1, t$V2, type="l",
         xlab="Distance (mm)",
         ylab="Absorbance",
         col="red",
         main=samples[i])
    abline(v=f$V1, col="lightgray", lty=3)
    text(f$V1, max(t$V2), labels=seq_along(f$V1), col="lightgray")
    dev.off()
}