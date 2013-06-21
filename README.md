#Bennett

+ high-vquest of joined reads
+ parse output for something meaningful

#Canine; Duval

+ go through email to determine comparison groups
+ script to run deseq across specified groups
+ attempt to characterize differentially expressed sequences

#Davidson; Human; T-Cell repertoire

#Deterding

#Hits-Clip

+ still need to run PK61-63
+ heatmap of miRNA, something that you can zoom into
+ matrix2png

#Poly(A)

+ SRCAP has a site that's 4 bases long, why?
+ get tracks up, change track heights, include shift directions
+ shortest height; session with all of the data; create groups
+ assign labels to tracks?
+ normals vs cancer; blah

+ intronic peaks
+ don't necessarily test with these sites, but may want to add yet another track
with all sites + intronic peaks
+ meme; motif finding; intronic memes for the clip data

+ may have to eventually compare present vs absent for peaks rather than relying
on dexseq calls for significance

+ lsf job submission help for austin
+ driver scripts

+ phyper
+ this set of miRNA in a given set of genes, what are the odds

+ call peaks using non-umi data
+ fake exons in order to visualize shifts

+ parse dexseq results into proximal, distal, no-change, no-test
+ merge bams, combine peaks, call classes, flag peaks not intersecting peaks as lower confidence

+ "in this test, this gene is observed to shift distally" -- using fold change
across significant sites of a single gene
+ figure out exactly how to characterize genes with more than 2 significant sites

+ map to mouse, then to human for peter's previous samples (samples with an 'a')
    + map some human to mouse and see how many map (percentage wise)
    + possible to mask human genome from mouse genome based on a certain length

    + testing mouse -- for each polya grab 25, 50, 75, 100 from the site
        + align to mouse and see what maps for each length
        + are there gains in using longer length reads

#Walter

+ levin classification; score
+ given list with direction, compute a score for each patient along with their group (TB, noinf, etc.)
+ calculate threshold for plot

+ adding total intensity at up-regulated transcripts and subtracting total intensity at all down-regulated transcripts
