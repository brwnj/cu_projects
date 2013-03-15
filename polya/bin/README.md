#poly(a)

pad the poly(a) sites by 10 bases. the poly(a) name field is poly(a)|gene
```
bedtools slop -l 5 -r 5 -i polya_3utr.bed -g hg18.sizes > polya_3_utr_slop.bed
```

the counts come from the bedgraphs
```
bedtools genomecov -strand + -bg -5 -ibam MP57.bam > MP57.pos.bedgraph
```

map the counts onto the poly(a) sites
```
bedtools map -c 4 -o sum -null 0 -a polya_3utr_slop.bed -b <sample>.bedgraph > <sample>.counts.bed
```

convert the counts to a friendlier format
```
awk 'BEGIN{OFS=FS="\t"}{split($4,name,"|"); print name[2],name[1],$7}' <sample>.counts.bed
```