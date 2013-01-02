.. _peterson:

******************************************************************************
Peterson
******************************************************************************

TODO
==============================================================================

* echo "gmap_build -g -d mm9 -k 15 chr*" | bsub-ez gmap


Reads
==============================================================================

Fastq QC

R1::

    PASS	Basic Statistics	        3Poly_1794e3_Mutant_NoIndex_L001_R1_001.fastq
    PASS	Per base sequence quality	3Poly_1794e3_Mutant_NoIndex_L001_R1_001.fastq
    PASS	Per sequence quality scores	3Poly_1794e3_Mutant_NoIndex_L001_R1_001.fastq
    WARN	Per base sequence content	3Poly_1794e3_Mutant_NoIndex_L001_R1_001.fastq
    WARN	Per base GC content	        3Poly_1794e3_Mutant_NoIndex_L001_R1_001.fastq
    FAIL	Per sequence GC content	    3Poly_1794e3_Mutant_NoIndex_L001_R1_001.fastq
    PASS	Per base N content	        3Poly_1794e3_Mutant_NoIndex_L001_R1_001.fastq
    PASS	Sequence Length Dist        3Poly_1794e3_Mutant_NoIndex_L001_R1_001.fastq
    PASS	Sequence Duplication Levels	3Poly_1794e3_Mutant_NoIndex_L001_R1_001.fastq
    PASS	Overrepresented sequences	3Poly_1794e3_Mutant_NoIndex_L001_R1_001.fastq
    WARN	Kmer Content	            3Poly_1794e3_Mutant_NoIndex_L001_R1_001.fastq

`Full QC Report R1`_

.. _Full QC Report R1: http://amc-einstein.ucdenver.pvt/~brownj/_static/3Poly_1794e3_Mutant_NoIndex_L001_R1_001_fastqc/fastqc_report.html

R2::
    
    PASS	Basic Statistics	        3Poly_1794e3_Mutant_NoIndex_L001_R2_001.fastq
    PASS	Per base sequence quality	3Poly_1794e3_Mutant_NoIndex_L001_R2_001.fastq
    PASS	Per sequence quality scores	3Poly_1794e3_Mutant_NoIndex_L001_R2_001.fastq
    WARN	Per base sequence content	3Poly_1794e3_Mutant_NoIndex_L001_R2_001.fastq
    WARN	Per base GC content	        3Poly_1794e3_Mutant_NoIndex_L001_R2_001.fastq
    FAIL	Per sequence GC content	    3Poly_1794e3_Mutant_NoIndex_L001_R2_001.fastq
    PASS	Per base N content	        3Poly_1794e3_Mutant_NoIndex_L001_R2_001.fastq
    PASS	Sequence Length Dist    	3Poly_1794e3_Mutant_NoIndex_L001_R2_001.fastq
    PASS	Sequence Duplication Levels	3Poly_1794e3_Mutant_NoIndex_L001_R2_001.fastq
    PASS	Overrepresented sequences	3Poly_1794e3_Mutant_NoIndex_L001_R2_001.fastq
    WARN	Kmer Content	            3Poly_1794e3_Mutant_NoIndex_L001_R2_001.fastq

`Full QC Report R2`_

.. _Full QC Report R2: http://amc-einstein.ucdenver.pvt/~brownj/_static/3Poly_1794e3_Mutant_NoIndex_L001_R2_001_fastqc/fastqc_report.html

Mapping
==============================================================================

``samtools flagstat`` output on GSNAP run:

===============     =========   =====
total reads         289312258
properly paired     186566127   64.4%
singletons          1496065     0.5%
===============     =========   =====

``bowtie2`` output::

    144656129 reads; of these:
      144656129 (100.00%) were paired; of these:
        8123049 (5.62%) aligned concordantly 0 times
        105831301 (73.16%) aligned concordantly exactly 1 time
        30701779 (21.22%) aligned concordantly >1 times
        ----
        8123049 pairs aligned concordantly 0 times; of these:
          1984713 (24.43%) aligned discordantly 1 time
        ----
        6138336 pairs aligned 0 times concordantly or discordantly; of these:
          12276672 mates make up the pairs; of these:
            7256924 (59.11%) aligned 0 times
            2566697 (20.91%) aligned exactly 1 time
            2453051 (19.98%) aligned >1 times
    97.49% overall alignment rate


Notes
==============================================================================

20120713
------------------------------------------------------------------------------

* project was delivered

20120601
------------------------------------------------------------------------------

* New project, detecting variant in chr4 region of mouse
* map with gsnap and bowtie and compare mapping stats
* first need to set up gsnap index