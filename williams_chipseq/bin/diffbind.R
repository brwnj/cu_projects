library(DiffBind)

setwd("~/projects/williams_chipseq/data/20140219")
williams = dba(sampleSheet="williams.csv")
williams = dba.count(williams)
williams = dba.contrast()