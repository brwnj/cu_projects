#library(ggplot2)

dt = read.table("~/projects/davidson/data/1/1.v_region_counts.txt", header=TRUE, sep="\t", row.names=1)
barplot(dt$coverage, names.arg=rownames(dt), las=3)