# Bennett
+ viewing in a tree may be helpful for manual picking
+ barcode to cell type conversion table!
+ markus really liked the local alignment results more so than the current way
+ figure out why CG1 is missing the match shown in 9/16 email
    + rerun alignments using local aligner and give those results
+ google doc to save numbers
+ from which population are the csf cells originating

#Cohrs
+ HOMER pipeline with design

#Duval

#HITS-CLIP
+ liftover to hg19
+ ucsc hub

#Leung
+ transfer all data
+ call variants on everything

# Poly(A)

### Mouse
+ map to mouse, then to human for peter's previous samples (samples with an 'a')
+ 86-...
+ map to human, map to mouse, what's the overlap
    + map some human to mouse and see how many map (percentage wise)
    + possible to mask human genome from mouse genome based on a certain length
    + testing mouse -- for each polya grab 25, 50, 75, 100 from the site
        + align to mouse and see what maps for each length
        + are there gains in using longer length reads


#Hanson
+ delivered

#Nicoli
+ liftover data to hg19
+ continue pipeline

#cpipe
+ filtering rRNA
+ we're seing 5' noise prior to the mIR sequence
+ the contaminant type will obviously depend on the library's composition
+ how do we know when the bin has crap on the edge of it
+ ends could be actually 3' or 5' real mutations
+ if something matches ribosomal RNA; continue
+ if illumina adapter is present, trim it
+ blat all sequences; prefiltering to rRNA database
+ align all sequences to genome, annotate to mapped position
+ annotate upstream and downstream to nearest gene as well
+ soft clipping will help trim those illumina adapters
+ vervet monkey just wants to know what the sequence is; annotation
+ crystal's project wants to know everything, including the garbage things
+ clean up counts by sno or mir or whatever
