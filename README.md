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

#Poly(A); hg18

+ call peaks, classify.
+ dexseq across samples using class 1 peaks inside of exons (including UTRs)
    + consensus of all peak regions at 10 bp spread from summit
+ classify proximal-distal using log2foldchange
+ group each classified dexseq results into categories

+ histone marks, ESTs to handle intergenic

+ cluster into patterns of states when comparing between replicates
+ all decreasing, all increasing, etc.

#Walter

+ scikit-learn for testing different classification methods
+ nick (k-nearest neighbors and random forests)
+ CART, svm, neural network
+ convert IDs from chauss... and finish implementing plotting along with corrplot-like plots
