#!/usr/bin/env python
# encoding: utf-8
"""
RNASEQ_PARAMETERS.py
"""
import cufflinks
import tophat
import os
import sys

##############################################################################
# DEFINE EXPERIMENT
##############################################################################
PAIRED_END = True
CONDITIONS = {'LT003990':{'read':'/mnt/storage3/brownj/lgrc/LT003990RU_R2_4_7-side_1.fastq /mnt/storage3/brownj/lgrc/LT003990RU_R2_4_7-side_2.fastq',
                          'case':'control'},
              'LT010012':{'read':'/mnt/storage3/brownj/lgrc/LT010012LU_R2_6_8-side_1.fastq /mnt/storage3/brownj/lgrc/LT010012LU_R2_6_8-side_2.fastq',
                          'case':'disease'}
             }
REFERENCE_SEQ = '/mnt/storage3/brownj/reference/Mus_musculus/Ensembl/NCBIM37/Sequence/WholeGenomeFasta/genome.fa'
GENE_ANNOTATION = '/mnt/storage3/brownj/reference/Mus_musculus/Ensembl/NCBIM37/Annotation/Genes/genes.gtf'

#QSUB = '/usr/local/bin/qsub'
#SERVER = 'genetics-compute.njrc.org'

# if reads are aligned, you may specify bams here
# list of bams as absolute paths
BAMS = []
# dictionary of transcripts full paths. key is label from CONDITIONS
TRANSCRIPTS = {}
MANIFEST = None

if len(BAMS) == 0:
    for sample, sample_info in CONDITIONS.iteritems():
        read = sample_info.get('read')
        if PAIRED_END:
            outputdir = '%s/%s' % (os.path.dirname(read.split(" ")[0]),
                    os.path.basename(read.split(" ")[0]).rstrip(".fastq"))
            bam = '%s/%s' % (outputdir, os.path.basename(read.split(" ")[0].rstrip(".fastq")))
            if os.path.exists(bam):
                BAMS.append(bam)
                sys.stderr.write("Bam found for: %s\n" % sample)
                continue
        else:
            outputdir = '%s/%s' % (os.path.dirname(read), 
                    os.path.basename(read).rstrip(".fastq"))
            bam = '%s/%s' % (outputdir, os.path.basename(read.rstrip(".fastq")))
            if os.path.exists(bam):
                BAMS.append(bam)
                sys.stderr.write("Bam found for: %s\n" % sample)
                continue
##############################################################################
# DEFINE TOPHAT OPTIONS
##############################################################################
        # Currently, the class with only read one option for PBS
        PBS_OPTIONS = ['#PBS -l walltime=72:00:00,nodes=1,ppn=8,mem=30gb']
        # Leave output dir as the first option.
        TOPHAT_OPTIONS = ['--output-dir %s' % outputdir,
                         '--mate-inner-dist 150',
                         '--num-threads 8',
                         '--segment-mismatches 2',
                         REFERENCE_SEQ.rstrip(".fa"),
                         read]
        th = RunTophat(PBS_OPTIONS, TOPHAT_OPTIONS)
        BAMS.append(th.runtophat())

# Ensure bams were added...
if len(BAMS) > 0:
    # ...then see if they're done
    for bam in BAMS:
        if os.path.exists(bam):
            label = ''
            for sample in CONDITIONS.keys():
                if sample in bam:
                    label = sample
            if label == '':
                sys.stderr.write("No label was assigned to %s. Skipping.\n" % bam)
                continue
##############################################################################
# DEFINE CUFFLINKS OPTIONS
##############################################################################
            # Currently, the class with only read one option for PBS
            PBS_OPTIONS = ['#PBS -l walltime=24:00:00,nodes=1,ppn=8,mem=30gb']
            # Leave output dir as the first option.
            CUFFLINKS_OPTIONS = ['--output-dir %s' % os.path.dirname(bam),
                                '--num-threads 8',
                                '--GTF %s' % GENE_ANNOTATION,
                                '--multi-read-correct',
                                '--label %s' % label,
                                '--frag-bias-correct %s' % REFERENCE_SEQ,
                                bam]
            cl = RunCufflinks(PBS_OPTIONS, CUFFLINKS_OPTIONS)
            TRANSCRIPTS[label] = cl.runcufflinks()
        else:
            sys.stderr.write("Bam %s doesn't exist. Is Tophat finished?\n" % bam)

if len(TRANSCRIPTS) > 0:
    transcriptpaths = []
    for sample, path in TRANSCRIPTS.iteritems():
        if not os.path.exists(path) and not os.path.getsize(path) > 1000:
            continue
        transcriptpaths.append(path)
    if len(transcriptpaths) > 0:
##############################################################################
# DEFINE CUFFMERGE OPTIONS
##############################################################################
        PBS_OPTIONS = ['#PBS -l walltime=48:00:00,nodes=1:ppn=8,mem=30gb']
        CUFFMERGE_OPTIONS = ['--ref-gtf %s' % GENE_ANNOTATION,
                            '--num-threads 8',
                            '--ref-sequence %s' % REFERENCE_SEQ,
                            manifest]
        cm = RunCuffmerge(transcriptpaths, PBS_OPTIONS, CUFFMERGE_OPTIONS)
        manifest = cm.createmanifest()
        cm.runcuffmerge()
CUFFMERGE = 'cuffmerge --ref-gtf %s --num-threads 8 --ref-sequence %s'
            % (GENE_ANNOTATION, REFERENCE_SEQ) #append the manifest.xml
CUFFDIFF = ''