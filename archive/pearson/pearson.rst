.. _pearson:

******************************************************************************
Pearson
******************************************************************************

QC
==============================================================================

* B1868_ATCACG_L001_R1_001_
* B1868_ATCACG_L001_R2_001_
* DisA1_TGACCA_L001_R1_001_
* DisA1_TGACCA_L001_R2_001_

.. _B1868_ATCACG_L001_R1_001: http://amc-sandbox.ucdenver.edu/~brownj/tesla/_static/pearson/B1868_ATCACG_L001_R1_001_fastqc/fastqc_report.html
.. _B1868_ATCACG_L001_R2_001: http://amc-sandbox.ucdenver.edu/~brownj/tesla/_static/pearson/B1868_ATCACG_L001_R2_001_fastqc/fastqc_report.html
.. _DisA1_TGACCA_L001_R1_001: http://amc-sandbox.ucdenver.edu/~brownj/tesla/_static/pearson/DisA1_TGACCA_L001_R1_001_fastqc/fastqc_report.html
.. _DisA1_TGACCA_L001_R2_001: http://amc-sandbox.ucdenver.edu/~brownj/tesla/_static/pearson/DisA1_TGACCA_L001_R2_001_fastqc/fastqc_report.html


Alignment
==============================================================================

gsnap aligner info::

    brownj at amc-tesla in ~/projects/pearson/results/20120816
    $ gsnap --version
    GSNAP version 2012-07-20 called with args: gsnap --version

    GSNAP: Genomic Short Nucleotide Alignment Program


Alignment stats
------------------------------------------------------------------------------

B1868

=============== =========== ===================    =============   ======================  ==========================
CATEGORY        TOTAL_READS PF_HQ_ALIGNED_READS    PF_INDEL_RATE   READS_ALIGNED_IN_PAIRS	PCT_READS_ALIGNED_IN_PAIRS
SECOND_OF_PAIR  50510392    49965174	            0.00035	        49986628	            0.989631 	
PAIR            101275898   100187808	            0.000331	    99973256	            0.987138
=============== =========== ===================    =============   ======================  ==========================

* not sure where the first of pair stats when. if you need them let me know...

DisA1

=============== =========== ===================  =============  ======================  ==========================
CATEGORY        TOTAL_READS PF_HQ_ALIGNED_READS  PF_INDEL_RATE  READS_ALIGNED_IN_PAIRS	PCT_READS_ALIGNED_IN_PAIRS
FIRST_OF_PAIR   34574198    33938980	         0.000618	    33398916	            0.966007	              
SECOND_OF_PAIR  35053427    34383546	         0.001025	    33398916	            0.9528	                  
PAIR            69627625    68322526	         0.000816	    66797832	            0.959358	              
=============== =========== ===================  =============  ======================  ==========================


20120815
==============================================================================
Candidate list::

   TTHERM_01308010 (Poc1)
   TTHERM_01084190 (Bbs1)
   TTHERM_00782070 (Bbs5)
   TTHERM_00161270 (Bbc31)
   TTHERM_00537420 (Fop1)
   TTHERM_00216010 (Ftt18)
   TTHERM_00068170 (Tcb2)


20120808
==============================================================================

* reference genome (tetrahymena) and gene annotation exist
* discover unique mutations of the mutant where wildtype will represent the background
* b = wildtype, control, subtract what's different from the reference
* d = mutant, pool of f2s
* data resides in: pilot/120727
* fastqc, align, gatk
* sb210 strain
* list of genes of interest
* macro vs mic genome reference mapping