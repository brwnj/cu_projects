# install.packages('corrplot')
# http://rstudio-pubs-static.s3.amazonaws.com/2107_4eb1adc1e4d44b93b6fde7eb801519fe.html

#normalized counts table
library(DESeq)
# python ../bin/data2matrix.py *MP*mirna_abundance.bed.gz > abundance_by_case.txt
ct = read.table("~/projects/hits-clip/data/abundance_by_case.txt", header=T, row.names=1)
cds = newCountDataSet(ct, conditions=c(rep("t",8),rep("c",8)))
cds = estimateSizeFactors(cds)
normcounts = counts(cds, normalized=T)

# correlation table
mat = cor(normcounts)
write.csv(mat, file="~/projects/hits-clip/data/correlation.csv")
write.csv((mat)^2, file="~/projects/hits-clip/data/correlation_coefficient.csv")

# plot
library(corrplot)
corrplot(mat, method="circle", order="alphabet", cl.lim=c(0,1))
corrRect(c(3,3,3,2,5), col="red", lwd=2)