#!/usr/bin/env python
# encoding: utf-8
"""
RunTophat.py

Created by Joe Brown on 2012-01-12.
"""

import sys
import os
import subprocess


class RunTophat(object):
    def __init__(self, pbsoptions, tophatoptions):
        super(RunTophat, self).__init__()
        self.pbsopts = pbsoptions
        self.tophatopts = tophatoptions


    def runtophat(self):
        outputdir = self.tophatopts[0].split(" ")[1]
        options = " ".join(o for o in self.tophatopts)
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        script = '%s/tophat.sh' % outputdir
        bam = '%s.bam' % self.tophatopts[-1].split(".")[0]
        sh = open(script, 'w')
        sh.write("#!/bin/sh\n")
        sh.write("%s\n" % self.pbsopts)
        sh.write("cd %s\n" % outputdir)
        sh.write("tophat %s\n" % options)
        sh.write("mv accepted_hits.bam %s" % bam)
        sh.close()
        subprocess.Popen(['qsub', script])
        # Needs to return the equivalent of an absolute path.
        return '%s/%s' % (outputdir, bam)