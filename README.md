#Bennett

+ take the most abundant sequence per UMI with highest quality for each direction
+ trim the primer sequences and annotate
+ join the reads
+ search for protein sequences

+ spreadsheet with protein sequence
+ annotate cdr based on protein sequence (allow up to 2 mismatches)
+ only interested in protein sequences containing YYC, YFC, and YSC

#Canine; Duval

TGGAATTCTCGGGTGCCAAGGAAC
TAGCTTATCAGACTGATGTTGACT
AAGCTGCCAGTTGAAGAACTGTTG

TGGAATTCTCGGGTGCCAAGGAAC
TAGCTTATCAGACTGATGTTGACT
TACCCTGTAGAACCGAATTTGTTGGAATTCTCGGG
AAGCTGCCAGTTGAAGAACTGTTGGAATTCTCGGG
TAGCTTATCAGACTGATGTTGATGGAATTCTCGGG
CAACGGAATCCCAAAAGCAGCTGTGGAATTCTCGG
TAGCAGCACGTAAATATTGGCGTGGAATTCTCGGG
TAGCTTATCAGACTGATGTTGACCTGGAATTCTCG
TGTCTGAGCGTCGCTTGGAATTCTCGGGTGCCAAG

#Davidson; Human; T-Cell repertoire

#Deterding

+ need to meet to discuss IPA

#Hits-Clip; hg18

+ still need to run PK61-63
+ heatmap of miRNA, something that you can zoom into
+ matrix2png

#Poly(A); hg18

+ add class to peak name

+ maybe just run dexseq without real replicates and see how many are significant
+ classify proximal-distal using log2foldchange
+ group each classified dexseq results into categories

+ map to mouse, then to human for peter's previous samples (samples with an 'a')
    + map some human to mouse and see how many map (percentage wise)
    + possible to mask human genome from mouse genome based on a certain length

    + testing mouse -- for each polya grab 25, 50, 75, 100 from the site
        + align to mouse and see what maps for each length
        + are there gains in using longer length reads

+ cluster into patterns of states when comparing between replicates
+ all decreasing, all increasing, etc.

#Walter

+ scikit-learn for testing different classification methods
+ nick (k-nearest neighbors and random forests)
+ CART, svm, neural network
+ convert IDs from chauss... and finish implementing plotting along with corrplot-like plots
