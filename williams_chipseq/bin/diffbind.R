setwd("~/Downloads/")
t = read.table("3B5_hela.AP2a.txt", sep="\t", header=TRUE)
plot(t[,1], t[,2], type="l", main="AP2a Positions in 3B5_hela", xlab="Distance From Peak (bp)", ylab="Sites per bp per Peak")
