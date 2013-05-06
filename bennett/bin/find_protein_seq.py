#!/usr/bin/env python
# encoding: utf-8
"""
7	LVESSGHY
9	RTLYGGIKDS 
14	GEPFSDSSGSYSLLHMQE
19	GLKSRNWYQGVDV
22	DRRTYYSGYQNLHFYGMDV
26	DPLSTYDSSGSYYDS
27	EGFISWEHYYSYHGMEV
32	DHSRIWWSEGWSRFSRGGQG
33	RSMVRGVNDV
89	RSMVRGVNDV
34	AMPTFGVAIIPGLDAGLSVIGVQV
36	AQEGSRGHALYVFHH
37	RTMGRGLNDY
18	RGLGRGVNDQ
41	RGLGRGVNDQ
43	GLRSSNWYEGMDV
44	SVQQVKVPTFDF
48	AATNWAGYYFDN
49	RSLVRGINDP
15	TNSLFGSSNRYYFNMDV
39	TNSLHGSAKIYLYSMDV
40	TNSLHGSAKIYLYSMDV
55	TNSLFGSSNRYFFNMDV
59	SNSLFGSPHRFYFNMDV
64	EIEF
65	EVIAMFGVVAGPAGGMDV
71	PGTDHGDFPCFDA
86	GHLDS
88	GLVFYDTLDS
91	QNLRLEQLHL
94	MNPDYGDFNYFDS
97	HLGPYLQWLLSGGPLSL

dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

In [57]: dna.translate()
Out[57]: Seq('MAIVMGR*KGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))
In [58]: dna[1:].translate()
Out[58]: Seq('WPL*WAAERVPD', HasStopCodon(ExtendedIUPACProtein(), '*'))
In [59]: dna[2:].translate()
Out[59]: Seq('GHCNGPLKGCPI', ExtendedIUPACProtein())
In [60]: dna[3:].translate()
Out[60]: Seq('AIVMGR*KGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))

"""
from Bio.Seq import Seq


def main(args):

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument()
    args = p.parse_args()
    main(args)
