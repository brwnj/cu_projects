#!/usr/bin/env python
# encoding: utf-8
"""
parse file names and move results to common. sample names must already exist
in common as directories.
"""
import os
import shutil

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def mv(src_dir, dest_dir):
    fileList = os.listdir(src_dir)
    for i in fileList:
        src = os.path.join(src_dir, i)
        dest = os.path.join(dest_dir, i)
        if os.path.exists(dest):
            if os.path.isdir(dest):
                mv(src, dest)
                continue
            else:
                os.remove(dest)
        shutil.move(src, dest_dir)


def main():
    # /vol1/home/brownj/projects/huang/results/common
    workdir = '/vol1/home/brownj/projects'
    commondir = 'results/common'
    # parse current working dir for project (always 5!)
    project = os.path.abspath(os.curdir).split("/")[5]
    # get sample names from dirs in common
    samples = os.listdir(os.path.join(workdir, project, commondir))
    # for each file in results, mv to sample dir
    for filename in os.listdir(os.curdir):
        for name in samples:
            if name in filename:
                # move the file from the current results dir to the common/sample dir
                shutil.copy(filename, os.path.join(workdir, project, commondir, name))


if __name__ == "__main__":
    main()
    
    
    
    



src_dir = '/Users/john.leschinski/Desktop/testSrc'
dest_dir = '/Users/john.leschinski/Desktop/testMove'
MoveOver(src_dir, dest_dir)