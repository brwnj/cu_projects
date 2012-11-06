#!/usr/bin/env python
# encoding: utf-8
"""
Created by Joe Brown on 2011-12-16.
"""
import os
import sys
import subprocess

class RunCuffmerge(object):
    """docstring for RunCuffmerge"""
    def __init__(self, transcripts, pbsoptions, cuffmergeoptions):
        super(RunCuffmerge, self, transcripts).__init__()
        self.transcripts = transcripts
        self.pbsopts = pbsoptions
        self.cuffmergeopts = cuffmergeoptions
        
    def createmanifest(self):
        """Writes the manifest of the transcripts.gtf files to include in
        Cuffmerge.
        """
        outputdir = os.path.dirname(os.path.dirname(self.transcripts[0]))
        manifest = '%s/manifest.txt' % outputdir
        sh = open(manifest, 'w')
        [sh.write('%s\n' % path) for path in self.transcripts]
        sh.close()
        return manifest


    def runcuffmerge(self):
        
                    sh = open(cuffmerge, 'w')
                    sh.write('#!/bin/sh\n')
                    sh.write('\n')
                    sh.write('cd %s\n' % outputdir)
                    sh.write('cuffmerge %s\n' % ())
                    sh.close()
                    Popen(['qsub', cuffmerge])
