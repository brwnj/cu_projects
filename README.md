#Bennett

+ spreadsheet with protein sequence
+ annotate cdr based on protein sequence (allow up to 2 mismatches)
+ only interested in protein sequences containing YYC, YFC, and YSC

#Davidson; Human; T-Cell repertoire

+ rerun assembly
+ rerun alpha
+ collapse beta, stats, etc...

#Deterding

+ gene list of the significant hits
+ run pathway analysis
+ blast hits across several popular species

#Hits-Clip; hg18

+ still need to run PK61-63
+ heatmap of miRNA, something that you can zoom into
matrix2png

#Poly(A); hg18

+ polyadb is 80% class 1, 20% class 3 when using the same read support cutoffs as the paper
+ 1440 (class 1 and 3) total peaks were classified using MP55
+ our class 1 consensus peaks (of 13 samples) peaks overlapped 75% of those 1440
    while containing another 8000 sites

+ dexseq across samples using class 1 peaks inside of exons (including UTRs)
+ classify proximal-distal using log2foldchange
+ group each classified dexseq results into categories

+ histone marks, ESTs to handle intergenic

+ map to mouse, then to human for peter's previous samples (samples with an 'a')
    + map some human to mouse and see how many map (percentage wise)
    + possible to mask human genome from mouse genome based on a certain length
    

+ cluster into patterns of states when comparing between replicates
+ all decreasing, all increasing, etc.

#Walter

+ scikit-learn for testing different classification methods
+ nick (k-nearest neighbors and random forests)
+ CART, svm, neural network
+ convert IDs from chauss... and finish implementing plotting along with corrplot-like plots
