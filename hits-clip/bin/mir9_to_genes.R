library(ggplot2)
library(reshape2)
library(plyr)

setwd("~/projects/hits-clip/data/20130328/")
df = read.table("GSE22219_highlow.txt", sep="\t", header=TRUE)
df <- subset(df, select = -c(genemap,mirsample,ER,Relapse,Time) )
# df = df[order(df$hsa.miR.9.5p),]
# mat = data.matrix(df)
# heatm = heatmap(mat, Rowv=NA, Colv=NA, col=cm.colors(256), scale="column")

df$Sample <- with(df, reorder(Sample, hsa.miR.9.5p))
dfm = melt(df)
dfm = ddply(dfm, .(variable), transform)
(p <- ggplot(dfm, aes(variable, Sample)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue"))
base_size <- 9
p + theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + opts(legend.position = "none", axis.ticks = theme_blank(), axis.text.x = theme_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))


vals = read.table("GSE22219_values.txt", sep="\t", header=TRUE)
vals = subset(vals, select = -c(genemap,mirsample,ER,Relapse,Time))
vals$Sample <- with(vals, reorder(Sample, hsa.miR.9.5p))
valsm = melt(vals)
valsm = ddply(valsm, .(variable), transform, rescale = rescale(value))
(p <- ggplot(valsm, aes(variable, Sample)) + geom_tile(aes(fill = rescale), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue"))
base_size <- 9
p + theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + opts(legend.position = "none", axis.ticks = theme_blank(), axis.text.x = theme_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))

# scaling fix
valss <- ddply(valsm, .(variable), transform, rescale = scale(value))
last_plot() %+% valss