# TODO

+ script to read bams, output normalized counts table
+ will require reference gtf, pybedtools, pandas, numpy

#Bennett

+ high-vquest of joined reads
+ parse output for something meaningful

#Canine; Duval

#Davidson; Human; T-Cell repertoire

#Deterding

#Hits-Clip

+ still need to run PK61-63
+ heatmap of miRNA, something that you can zoom into
+ matrix2png

#Poly(A)

+ get tracks up, change track heights, include shift directions
+ shortest height; session with all of the data; create groups
+ assign labels to tracks?
+ normals v cancer

+ intronic peaks; start and stop of gene model as separate run; let austin know its location
+ don't necessarily test with these sites, but may want to add yet another track with all sites + intronic peaks

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
