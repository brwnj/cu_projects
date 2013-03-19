#!/usr/bin/env python
# encoding: utf-8
"""
reverse complement
"""

def reverse(s):
    """return in reverse order"""
    seq = list(s)
    seq.reverse()
    return ''.join(seq)
 
def complement(s):
    """return complementary seq"""
    basecomplement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N',
    'a':'t', 'c':'g', 'g':'c', 't':'a'}
    seq = list(s)
    seq = [basecomplement[b] for b in seq]
    return ''.join(seq)
 
def reversecomplement(s):
    """return reverse complement"""
    s = reverse(s)
    s = complement(s)
    return s

gencode = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
	'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

def valid_umi(iupac, umi):
    """parse UMI sequence to validate against IUPAC sequence."""
    IUPAC_definitions = {"A":"A","T":"T","C":"C","G":"G","R":"GA","Y":"TC",
                            "M":"AC","K":"GT","S":"GC","W":"AT","H":"ACT",
                            "B":"GTC","V":"GCA","D":"GAT","N":"GATC"}
    for code, base in izip(iupac, umi):
        try:
            if not base in IUPAC_definitions[code]:
                return False
        except KeyError:
            return False
    return True

def main(args):
    print reversecomplement(args.seq)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("seq")
    main(p.parse_args())

# @M00658:32:000000000-A3ARU:1:1101:14184:1325 1:N:0:9
# NGGTGCCAGGGGGAAGACCGATGGGCCCTTGGTGGAGGCTGATGAGACGGTGACCAGGGTTCCCTGGCTCCAGTAGTCAAAGAAGTACCAGCCACTGCTCTGGTGCGGATCTTTCGCACAGTAATATACGGCCGTGTCCTCGGCTCTCAGGCTGTTCACTTGCAGATACAGCGTGCTCTTGGAATTGTCTCTGGAGATGGTGAACCGGCCCTTCACGGAGTCTGCGTAGTATGTGCTACCACCACTAGCAC
# +
# #>>AA>FBCAADEGGGGGGEAEEECFGGHHHFDFFHHFFHG1FFGCEGEGGHFHHHHGGHHHHGHFAGGGGFFHHHHHGHHHHHHGHHHHHHGHHHHHHHHHFGHGGGGGGHHGHGGGGGHHHHFHHHHGGGGGC?GHFDGGDEFGHHHFFFHHHHHHHHGHGFHEFCGHHHGDADGGEBFFFFFGF;;9CFBBFFFGFBCFFGB.-9@->./.;9A.>DDFFFFF..BAFFFFFFBFFFBF.;?FFBFB/
# 
# @M00658:32:000000000-A3ARU:1:1101:14184:1325 2:N:0:9
# CAGTGTCAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTTTGCCATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCTTTTATTAGTGCTAGTGGGGGGAGCACATACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAAAGCACGCTGTATCTGCA
# +
# BBBB@FFFFFFFBGGFGGG4GGHHFHGFHB?EEGGEHGGGFFFFE3FFFHGCEFGBC02G00DFDGDGBGHHHGHHHGF01GGH1GBGHHHHHFHGFHHHGH0DHHHFHHFGHHHGHF?@GFGG?..AFBFDE.AEDGA-ED.:9FFFA.99FFFF/B//////9/B//;.;B---9...9F//9:ABD--9;FBA.;.EFFF.;9B-..9F/9F//9////.9.BFFFF/B//.9.;.--9:ABB//9/F