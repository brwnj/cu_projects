library('DiffBind')
setwd(system.file("extra", package="DiffBind"))
data(tamoxifen_peaks)
plot(tamoxifen)
data(tamoxifen_counts)
plot(tamoxifen)
data(tamoxifen_analysis)
plot(tamoxifen)
dba.plotMA(tamoxifen)
dba.plotMA(tamoxifen, bXY=T)
#looking at all binding sites
dba.plotPCA(tamoxifen, contrast=1,th=1)
#only differentially bound sites, fdr<0.1
dba.plotPCA(tamoxifen,contrast=1, th=.1)
#pca based on tissue type specified in meta
dba.plotPCA(tamoxifen, DBA_TISSUE)
pvals = dba.plotBox(tamoxifen)
corvals = dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE)
#dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE, scale=row)

data(tamoxifen_peaks)
olap.rate = dba.overlap(tamoxifen,mode=DBA_OLAP_RATE)
plot(olap.rate,type='b')
#peaks found in each of these biological replicate
dba.plotVenn(tamoxifen,tamoxifen$masks$MCF7)

#choosing an overlap rate per group
dba.overlap(tamoxifen,tamoxifen$masks$Resistant,mode=DBA_OLAP_RATE)
dba.overlap(tamoxifen,tamoxifen$masks$Responsive,mode=DBA_OLAP_RATE)
tamoxifen = dba.peakset(tamoxifen,tamoxifen$masks$Resistant, sampID="Resistant",minOverlap=2)
tamoxifen = dba.peakset(tamoxifen,tamoxifen$mask$Responsive, sampID="Responsive",minOverlap=3)
dba.plotVenn(tamoxifen,tamoxifen$masks$Consensus)
