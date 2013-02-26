#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division 
from os.path import splitext

from cogent.util.option_parsing import parse_command_line_parameters, make_option
from cogent import LoadSeqs, DNA
from cogent.parse.fasta import MinimalFastaParser
from qiime.util import FunctionWithParams

import cogent.app.muscle_v38
import cogent.app.clustalw
import cogent.app.mafft
import cogent.app.raxml_v730
import cogent.app.fasttree
import cogent.app.fasttree_v1
import cogent.app.clearcut
import warnings

class TreeBuilder(FunctionWithParams):
    """A TreeBuilder takes a aligned set of sequences and returns a tree.

    This is an abstract class: subclasses should implement the __call__
    method.

    Note: sequence ids should be preserved during this process, i.e. the
    description lines should be saved/restored if the phylogeny app is
    destructive to them.

    Specific wrappers need to be written for nontraditional approaches,
    including:

    (a) RAxML assignment, where sequences are assigned to internal nodes of
        an existing tree
    (b) BLAST-based assignment, where sequences are assigned to existing
        nodes of a tree based on best blast hit (assigned only to terminal
        nodes if using a single hit, but could easily imagine assigning to
        an internal node based on a set of indistinguishably good hits).
    (c) BLAST-like approach using VMATCH or other suffix array library, or 
        using oligonucleotide freqs like the RDP classifier does to assign
        to an arbitrary internal node in an existing tree, etc.
    """
    Name = 'TreeBuilder'

    def __init__(self, params):
        """Return new TreeBuilder object with specified params.
        
        Note: expect params to contain both generic and per-method (e.g. for
        raxml vs. fasttree vs. whatever) params, so leaving it as a dict 
        rather than setting attributes. Some standard entries in params are:

        Application: 3rd-party application used, if any, e.g. raxml
        """
        self.Params = params

    def __call__ (self, aln_path, result_path=None, log_path=None):
        """Returns tree from alignment.
        
        Parameters:
        aln_path: path to file of aligned sequences
        result_path: path to file of results. If specified, should
        dump the result to the desired path as fasta, otherwise should
        return cogent.core.alignment.DenseAlignment object.
        log_path: path to log, which should include dump of params.
        """
        raise NotImplementedError, "TreeBuilder is an abstract class"

class CogentTreeBuilder(TreeBuilder):
    """Generic tree builder using Cogent tree methods."""

    Name = 'CogentTreeBuilder'

    def getResult(self, aln_path, *args, **kwargs):
        """Returns alignment from sequences.
        
        Currently does not allow parameter tuning of program and uses
        default parameters -- this is bad and should be fixed.

        #TODO: allow command-line access to important aln params.
        """
        module = self.Params['Module']
        # standard qiime says we just consider the first word as the unique ID
        # the rest of the defline of the fasta alignment often doesn't match
        # the otu names in the otu table
        seqs = LoadSeqs(aln_path, Aligned=True, label_to_name=lambda x: x.split()[0])
        result = module.build_tree_from_alignment(seqs, moltype=DNA)

        try: 
            root_method = kwargs['root_method']
            if root_method == 'midpoint':
                result = root_midpt(result)
            elif root_method == 'tree_method_default':
                pass
        except KeyError:
            pass
        return result

    def __call__(self, result_path=None, log_path=None, *args, **kwargs):
        """Calls superclass method to align seqs"""
        return FunctionWithParams.__call__(self, result_path=result_path,
                                           log_path=log_path, *args, **kwargs)


tree_method_constructors = {}
tree_module_names = {'muscle':cogent.app.muscle_v38,
    'clustalw':cogent.app.clustalw,
    #'mafft':cogent.app.mafft,   
    # current version of Mafft does not support tree building
    'fasttree':cogent.app.fasttree,'fasttree_v1':cogent.app.fasttree_v1,
    'raxml_v730':cogent.app.raxml_v730,'clearcut':cogent.app.clearcut
    }

def root_midpt(tree):
    """ this was instead of PhyloNode.rootAtMidpoint(), which is slow and broke
    
    this should be deprecated in a future release once the release version
    of PyCogent's tree.rootAtMidpoint() is identical to this function
    
    this fn doesn't preserve the internal node naming or structure,
    but does keep tip to tip distances correct.  uses unrootedDeepcopy()
    """
    #max_dist, tip_names, int_node = getMaxTipTipDistance(tree)
    max_dist, tip_names, int_node = tree.getMaxTipTipDistance()

    half_max_dist = max_dist/2.0
    if max_dist == 0.0: # only pathological cases with no lengths
        return tree.unrootedDeepcopy()
    tip1 = tree.getNodeMatchingName(tip_names[0])
    tip2 = tree.getNodeMatchingName(tip_names[1])
    lca = tree.getConnectingNode(tip_names[0],tip_names[1]) # last comm ancestor
    if tip1.distance(lca) > half_max_dist:
        climb_node = tip1
    else:
        climb_node = tip2
        
    dist_climbed = 0.0
    while dist_climbed + climb_node.Length < half_max_dist:
        dist_climbed += climb_node.Length
        climb_node = climb_node.Parent
    
    # now midpt is either at on the branch to climb_node's  parent
    # or midpt is at climb_node's parent
    # print dist_climbed, half_max_dist, 'dists cl hamax'
    if dist_climbed + climb_node.Length == half_max_dist:
        # climb to midpoint spot
        climb_node = climb_node.Parent
        if climb_node.isTip():
            raise RuntimeError('error trying to root tree at tip')
        else:
            # print climb_node.Name, 'clmb node'
            return climb_node.unrootedDeepcopy()
        
    else:
        # make a new node on climb_node's branch to its parent
        tmp_node_name = "TMP_ROOT_NODE_NAME"
        parent = climb_node.Parent
        parent.removeNode(climb_node)
        climb_node.Parent = None
        new_node = parent.__class__()
        new_node.Name = tmp_node_name

        # adjust branch lengths
        old_br_len = climb_node.Length
        climb_node.Length = half_max_dist - dist_climbed
        new_node.Length = old_br_len - climb_node.Length

        if climb_node.Length < 0.0 or new_node.Length < 0.0:
            raise RuntimeError('attempting to create a negative branch length!')

        # re-attach tree
        parent.append(new_node)
        new_node.append(climb_node)
        
        # reroot and remove the temporary node name
        new_tree = tree.rootedAt(tmp_node_name)
        new_root = new_tree.getNodeMatchingName(tmp_node_name)
        new_root.Name = None
        
        return new_tree


script_info={}
script_info['brief_description']="""Make Phylogeny"""
script_info['script_description']="""Many downstream analyses require that the phylogenetic tree relating the OTUs in a study be present. The script make_phylogeny.py produces this tree from a multiple sequence alignment. Trees are constructed with a set of sequences representative of the OTUs, by default using FastTree (Price, Dehal, & Arkin, 2009)."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""A simple example of make_phylogeny.py is shown by the following command, where we use the default tree building method (fasttree) and write the file to the current working directory without a log file:""","""%prog -i $PWD/aligned.fasta -o $PWD/rep_phylo.tre"""))
script_info['script_usage'].append(("""""","""Alternatively, if the user would prefer using another tree building method (i.e. clearcut (Sheneman, Evans, & Foster, 2006)), then they could use the following command:""","""%prog -i $PWD/aligned.fasta -t clearcut"""))
script_info['output_description']="""The result of make_phylogeny.py consists of a newick formatted tree file (.tre) and optionally a log file. The tree file is formatted using the Newick format and this file can be viewed using most tree visualization tools, such as TopiaryTool, FigTree, etc.

The tips of the tree are the first word from the input sequences from the fasta file, e.g.: '>101 PC.481_71 RC:1..220' is represented in the tree as '101'."""
script_info['required_options']=[
    make_option('-i','--input_fp',action='store',
     type='existing_filepath',dest='input_fp',help='Path to read '+\
     'input fasta alignment, only first word in defline will be considered')
]
valid_root_methods = ['midpoint','tree_method_default']

script_info['optional_options']=[\
    make_option('-t','--tree_method',action='store',type='choice', choices=list(tree_module_names.keys()),
          help='Method for tree building. Valid choices are: '+\
          ', '.join(tree_module_names.keys())+\
          ' [default: %default]', default='fasttree'),
          
    make_option('-o','--result_fp',action='store',type='new_filepath',
          help='Path to store '+\
          'result file [default: <input_sequences_filename>.tre]'),
          
    make_option('-l','--log_fp',action='store',type='new_filepath',
          help='Path to store '+\
          'log file [default: No log file created.]'),
          
    make_option('-r','--root_method',action='store',type='choice', choices=list(valid_root_methods),
        help='method for choosing root of phylo tree'+\
        '  Valid choices are: '+ ', '.join(valid_root_methods) +\
        ' [default: tree_method_default]',
        default='tree_method_default'),
]
script_info['version'] = "1.0"

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

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
    if not result_path: # empty or None
        fpath, ext = splitext(input_seqs_filepath) # fpath omits extension
        result_path = fpath + ".tre"
    
    open(result_path,'w').close() # touch
    log_path = opts.log_fp
    if log_path != None:
        open(log_path,'w').close()

    if tree_builder_type=='Constructor':
        tree_builder(input_seqs_filepath,\
        result_path=result_path,log_path=log_path,failure_path=failure_path)
    elif tree_builder_type=='Cogent':
        tree_builder(result_path, aln_path=input_seqs_filepath,
            log_path=log_path, root_method=opts.root_method)


if __name__ == "__main__":
    main()