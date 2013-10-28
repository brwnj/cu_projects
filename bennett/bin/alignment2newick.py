#!/usr/bin/env python
# encoding: utf-8
"""
Convert MUSCLE alignment fasta to Newick string (.tre)
"""
import argparse

from qiime.make_phylogeny import tree_module_names, tree_method_constructors, CogentTreeBuilder

def main(fasta, tre):
    if not (opts.tree_method in tree_method_constructors or
            opts.tree_method in tree_module_names):
        option_parser.error(\
         'Invalid alignment method: %s.\nValid choices are: %s'\
         % (opts.tree_method,\
            ' '.join(tree_method_constructors.keys() +
                tree_module_names.keys())))
    try:
        tree_builder_constructor =\
            tree_method_constructors[opts.tree_method]
        tree_builder_type = 'Constructor'
        params = {}
        tree_builder = tree_builder_constructor(params)
    except KeyError:
        tree_builder = CogentTreeBuilder({
                'Module':tree_module_names[opts.tree_method],
                'Method':opts.tree_method
                })
        tree_builder_type = 'Cogent'

    input_seqs_filepath = opts.input_fp
    result_path = opts.result_fp
    
    open(result_path,'w').close()
    tree_builder(result_path, aln_path=input_seqs_filepath, log_path=log_path, root_method=opts.root_method)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("fasta", help="aligned sequences")
    p.add_argument("tre", help="tre file name")
    args = p.parse_args()
    main(args.fasta, args.tre)
