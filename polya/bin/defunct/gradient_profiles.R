library(ggplot2)
setwd("~/projects/polya/data/02052013/")

df = read.table("new_parse_method_1.tab", header=TRUE, sep="\t")

need a rectangle from 0 to the last distance of fraction == 0

rect_left = c()
rect_right = c()

positions = c(1,3,5,7,9,11,13,15,17,19,21)
for(p in seq_along(positions)){
    print(p)
    print(min(df[df$Fraction %in% positions[p], "Distance"]))
    print(max(df[df$Fraction %in% positions[p], "Distance"]))
}

min(df[df$Fraction %in% positions[5], "Distance"])



p <- ggplot(df, aes(x=Distance, y=Absorbance)) + 
    geom_line() +
    theme_bw()
p

p



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