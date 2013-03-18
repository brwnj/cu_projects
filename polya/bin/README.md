#poly(a)

pad the poly(a) sites by 10 bases. the poly(a) name field is poly(a)|gene
```
bedtools slop -l 5 -r 5 -i polya_3utr.bed -g hg18.sizes > polya_3_utr_slop.bed
```