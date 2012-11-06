#!/usr/bin/env python
# encoding: utf-8
"""
Runs Cufflinks on supplied mapped reads.

Created by Joe Brown on 2011-12-12.
"""
import os
import sys
import subprocess


class RunCufflinks(object):
    """Runs cufflinks as per options defined in user parameters file."""
    def __init__(self, pbsoptions, cufflinksoptions):
        super(RunCufflinks, self, pbsoptions, cufflinksoptions).__init__()
        self.pbsopts = pbsoptions
        self.cufflinksopts = cufflinksoptions


    def runcufflinks(self):
        outputdir = self.cufflinksopts[0].split(" ")[1]
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        script = '%s/cufflinks.sh' % outputdir
        sh = open(script, 'w')
        sh.write('#!/bin/sh\n')
        sh.write('%s\n' % self.pbsopts)
        sh.write("cd %s\n" % outputdir)
        sh.write("cufflinks %s\n" % options)
        sh.close()
        subprocess.Popen(['qsub', script])
        return '%s/transcripts.gtf' % outputdir