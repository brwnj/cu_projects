#Bennett

+ You want numbers of unique protein sequences per barcode for patient 1
```
for f in *.aa.txt; do echo $f; awk '{if(NR==1){print};if($3=="productive"){split($2,n,":"); print n[2]":"n[3]":"$15,$0}}' $f | sort -u -k1,1 | wc -l; done
```
```
for f in *.aa.txt; do awk '$3=="productive"{split($2,n,":"); print n[2]":"n[3]":"$15,$0}' $f | sort -u -k1,1 | cut -f2- > ${f%.aa*}.unique.headless.txt; done
```
+ per umi filtering...

+ take the most abundant sequence per UMI with highest quality for each direction
+ trim the primer sequences and annotate
+ join the reads
+ search for protein sequences

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
