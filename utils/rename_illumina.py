#! /usr/bin/env python
# encoding: utf-8
"""
rename fastqs named via illumina naming convention. assumes everything after 
first period should remain.
"""
import os

def main(args):
    for file in args.files:
        full_name, ext = file.split(".", 1)
        sample_name = "_".join(full_name.split(args.delimiter)[:args.save])
        new_name = "%s.%s" % (sample_name, ext)
        os.rename(file, new_name)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+", help="files to rename")
    req = p.add_argument_group("required arguments")
    req.add_argument("-s", "--save", type=int, required=True, 
            help="number of name pieces to save after splitting by delimiter")
    p.add_argument("-d", "--delimiter", default="_", 
            help="delimiter used in file name [ _ ]")
    main(p.parse_args())