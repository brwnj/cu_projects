#!/usr/bin/env python
# encoding: utf-8
"""
fqreader.py
"""

import sys
import os
import unittest
from FastQRecord import FastQRecord


class FastQReader:
    def __init__(self, fastq_file):
        self.fastq = fastq_file

    def __iter__(self):
        return self

    def next(self):
        name = self.fastq.next()
        sequence = self.fastq.next()
        comment = self.fastq.next()
        quality = self.fastq.next()
        if not quality.endswith("\n"):
            quality += "\n"
        if not name.startswith("@"):
            raise ValueError("ID line does not start with '@'. Out of sync or invalid file type.")
        if not comment.startswith("+"):
            raise ValueError("Comment line does not start with '+'. Out of sync or invalid file type." )
        if len(comment) > 2 and name[1:] != comment[1:]:
            raise ValueError( "Name and comment lines are not equal." )
        return FastQRecord(name[1:-1], sequence[:-1], comment[1:-1], quality[:-1])


class FastQReaderTests(unittest.TestCase):
    def setUp(self):
        pass


if __name__ == '__main__':
    unittest.main()