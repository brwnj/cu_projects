#! /usr/bin/env python
# encoding: utf-8
"""
Rename fastqs named via illumina naming convention. File name is split by
delimiter. Sections to save are numbered as 1-based.
"""
import os

def main(args):
    rename = {}
    new_names = []
    for file in args.files:
        full_name, ext = file.split(".", 1)
        name_pieces = full_name.split(args.delimiter)
        save = args.save.split(",")
        sample_name = "_".join([name_pieces[int(s)-1] for s in save])
        new_name = "%s.%s" % (sample_name, ext)
        new_names.append(new_name)
        rename[file] = new_name
    new_names = set(new_names)
    # catch non-unique names
    assert(len(new_names) == len(args.files))
    for file, new_name in rename.iteritems():
        os.rename(file, new_name)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+", help="files to rename")
    req = p.add_argument_group("required arguments")
    req.add_argument("-s", "--save", required=True, 
            help="pieces of name to save by comma separated, 1-based number list")
    p.add_argument("-d", "--delimiter", default="_", 
            help="delimiter used in file name [ _ ]")
    main(p.parse_args())