#!/usr/bin/env python
# encoding: utf-8
"""
fqrecord.py
"""

import sys
import os
import unittest


class FastQRecord:
    def __init__(self, name, sequence, comment, quality):
        self.name = name
        self.sequence = sequence
        self.comment = comment
        self.quality = quality


class FastQRecordTests(unittest.TestCase):
    def setUp(self):
        pass


if __name__ == '__main__':
    unittest.main()