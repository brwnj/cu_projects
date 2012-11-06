library(seqbias)
library(Rsamtools)
library(ggplot2)
ref_fn <- system.file( "extra/example.fa", package = "seqbias" )
ref_f <- FaFile( ref_fn )
open.FaFile( ref_f )
reads_fn <- system.file( "extra/example.bam", package = "seqbias" )

#sampling
ref_seqs <- scanFaIndex( ref_f )
I <- random.intervals( ref_seqs, n = 5, m = 100000 )
#sequences
seqs <- scanFa( ref_f, I )
neg_idx <- as.logical( I@strand == '-' )
seqs[neg_idx] <- reverseComplement( seqs[neg_idx] )
#counts
counts <- count.reads( reads_fn, I, binary = T )
#frequencies
freqs <- kmer.freq( seqs, counts )
if( require(ggplot2) ) {
    P <- qplot(x=pos,y = freq,ylim = c(0.15,0.4), color = seq, data = freqs, geom = "line" )
    P <- P + facet_grid( seq ~ . )
    print(P)
 }
else {
   par( mar = c(5,1,1,1), mfrow = c(4,1) )
   with( subset( freqs, seq == "a" ),
 plot( freq ~ pos, ylim = c(0.15,0.4), sub = "a", type = 'l' ) )
   with( subset( freqs, seq == "c" ),
 plot( freq ~ pos, ylim = c(0.15,0.4), sub = "c", type = 'l' ) )
   with( subset( freqs, seq == "g" ),
 plot( freq ~ pos, ylim = c(0.15,0.4), sub = "g", type = 'l' ) )
   with( subset( freqs, seq == "t" ),
 plot( freq ~ pos, ylim = c(0.15,0.4), sub = "t", type = 'l' ) )
}
#compensation
#training
sb <- seqbias.fit( ref_fn, reads_fn, L = 25, R = 45 )
#prediction
bias <- seqbias.predict( sb, I )
#adjustment
counts.adj <- mapply( FUN = `/`, counts, bias, SIMPLIFY = F )
freqs.adj <- kmer.freq( seqs, counts.adj )
par( mar = c(5,1,1,1), mfrow = c(4,1) )
with( subset( freqs.adj, seq == "a" ),
 plot( freq ~ pos, ylim = c(0.15,0.4), sub = "a", type = 'l' ) )
   with( subset( freqs.adj, seq == "c" ),
 plot( freq ~ pos, ylim = c(0.15,0.4), sub = "c", type = 'l' ) )
   with( subset( freqs.adj, seq == "g" ),
 plot( freq ~ pos, ylim = c(0.15,0.4), sub = "g", type = 'l' ) )
   with( subset( freqs.adj, seq == "t" ),
 plot( freq ~ pos, ylim = c(0.15,0.4), sub = "t", type = 'l' ) )

