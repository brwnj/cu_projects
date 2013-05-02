#Davidson; Human; T-Cell repertoire

+ rerun assembly
+ rerun alpha
+ collapse beta, stats, etc...

#Deterding

+ gene list of the significant hits
+ run pathway analysis
+ blast hits across several popular species

##What has been done
+ pair down number of proteins (1200 -> 100 disease specific)
+ panther and reactome for pathway analysis; access to gene-go
+ connectivity maps for drug interaction
+ down regulation is almost entirely novel
+ within CF, can you classify groups?
+ are the groups clinically different; they may have different bug
+ which proteins are responsible for the differences
+ moderate (35) vs severe (14) disease; total n = 
+ phenotypic data is abundant
+ FEV measurement to determine classifying protein for rapid declining FEV

#Hits-Clip; hg18

+ still need to run PK61-63

#Poly(A); hg18

+ call peaks, classify.
+ dexseq across samples using class 1 peaks inside of exons (including UTRs)
+ classify proximal-distal using log2foldchange
+ group each classified dexseq results into categories

+ histone marks, ESTs to handle intergenic

+ dexseq across samples, ctrl, estr tx, tam
+ cluster into patterns of states when comparing between replicates
+ all decreasing, all increasing, 

+ something about noise of mouse cells polluting alignments

#Walter

+ scikit-learn for testing different classification methods
+ nick (k-nearest neighbors and random forests)
+ CART, svm, neural network
+ convert IDs from chauss... and finish implementing plotting along with corrplot-like plots
