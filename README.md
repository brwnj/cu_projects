#Bennett

+ take the most abundant sequence per UMI with highest quality for each direction
+ trim the primer sequences and annotate
+ join the reads
+ search for protein sequences

+ spreadsheet with protein sequence
+ annotate cdr based on protein sequence (allow up to 2 mismatches)
+ only interested in protein sequences containing YYC, YFC, and YSC

#Canine; Duval

#Davidson; Human; T-Cell repertoire

#Deterding

#Hits-Clip

+ still need to run PK61-63
+ heatmap of miRNA, something that you can zoom into
+ matrix2png

#Poly(A)

+ add class to peak name
+ update ucsc
+ run dexseq for all combinations (1v1 with fake rep)
+ send dexseq results

+ start work on fisher testing of pairs to see whether we're getting the right answer

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
+ convert IDs from chauss... and finish implementing plotting along with corrplot-like plots

+ given list with direction, compute a score for each patient along with their group (TB, noinf, etc.)
+ calculate threshold for plot
