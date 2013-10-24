#!/usr/bin/env python
# encoding: utf-8

from editdist import distance
from itertools import combinations
from collections import defaultdict

seqs = "AEETRFTITTSFDP AGETRFTVTTSFDP AGGTRFTVTTSFDP AIETMFTVTTSVDP AKEPRFTVPTSFDP AKEPRFTVTTSFDP AKETRFTVTTSFDP ARATRVTVTTPFDP ARDIRFTFTTSIDP ARDTRFTVTTSVDP AREARFTVTTSFDP AREPRFTVTTSFDP ARESRFTVTTSFDP ARETGFTVTTSFDP ARETKFTVTTSFDP ARETQFTVTTSFHP ARETRCTVTTSFDP ARETRFSVTTSFDP ARETRFTVATSFDP ARETRFTVPTSFDP ARETRFTVSTSFDP ARETRFTVSTSLDA ARETRFTVTNSFDP ARETRFTVTSSFDA ARETRFTVTSSFDP ARETRFTVTTPFDP ARETRFTVTTSCDL ARETRFTVTTSFDL ARETRFTVTTSFDP ARETRFTVTTSFGP ARETRFTVTTSIDP ARETRLTVTTSFDP ARETRSTVTTSFDP ARETRVTVTTSFDP ARGTRFTVTTSFDP ARQTRFTVTTSFDP ERETRFTVTTSFDP SRETRFTVTTSFDP TRETRFTVTTSFDP".split(" ")
# distances = defaultdict(list)
# for (a,b) in combinations(seqs, 2):
#     d = distance(a,b)
#     distances[a].append(d)
#     distances[b].append(d)
# # gives nice lists that allow you to quickly find which one has the most 1s
# # but it isn't very nice for collapsing
# print distances
print "\n".join(seqs)