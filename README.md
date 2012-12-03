==============================================================================
Peaktools
==============================================================================

prints None sometimes:
    
    chrY	None	None	MDX_22_AGTTCC_L003_R1_001_pos	None	+

==============================================================================
Artinger; Zebrafish; ChIP-Seq
==============================================================================

==============================================================================
Hits-Clip
==============================================================================

* Align with bowtie suppressing anything that maps to more than one place
* Filter the bam of anything that aligned overlapping the regions of known rRNAs
* Deliver tracks

* ayb with shorter and shorter (less priority)

==============================================================================
Jacobsen; Human; ChIP-Seq
==============================================================================

To elucidate endogenous binding sites for progesterone receptor-A (PRA) we used
three different breast cancer cell sublines of T47D cells; PR-negative cells
(Y cells), cells that constitutively express PRA (YA cells) and cells that 
inducibly express PRA upon treatment with ponA (YiA cells). These cells were 
treated with vehicle (EtOH), Progesterone (P4), or an anti-progestin (ZK). We 
wish to examine the binding of PR to DNA in the absence of any ligand (EtOH), 
the presence of progesterone (P4), or the presence of anti-progestin (ZK). 
Labels for the samples are as follows: Y+EtOH YA+EtOH YA+P4 YA+ZK
(1 hour treatment) YA+ZK (24 hour treatment) YiA+EtOH YiA+ponA+EtOH 
(cells induced to express PRA) YiA+ponA+ZK (1 hour) YiA+ponA+ZK 
(24 hour treatment) YiA+ponA+P4 While the cell sublines are different, the 
two PRA+ cells should be similar with regard to PR binding to target 
sequences. For example, YA+EtOH should be similar to YiA+ponA+EtOH, YA+P4 
should be similar to YiA+ponA+P4, YA+ZK 1h is similar to YiA+PonA+ZK 1h, 
YiA EtOH is similar to Y EtOH etc. One caveat is that because PRA is inducibly 
expressed in YiA cells, the cells make "leak" and express some PRA 
(even in the sample only treated with EtOH). Is it possible to analyze the 
samples first as replicates then in singlicate?

==============================================================================
Poly(A)
==============================================================================

    samtools view bam | python dexseq_count.py flattened.gff - out.counts
    rm *x.counts
    for f in *.counts; do awk 'BEGIN{OFS=FS="\t"}!/^_/{foo=int(rand()*10);if($2!=0){print $1,$2+foo}else{print}}' $f > $(basename $f .counts)x.counts; done

* neg and pos split everything; annotation, bam, etc
* tracks have to be shown only using 5' read ends
* re-run dexseq

==============================================================================
Walter
==============================================================================

* rerun analysis on human samples
* find where some of the sequences are that are DE in H37Rv -- they be DE in hg19
* annotate to genome feature, primary interest is 3' UTR