#Bennett

+ You want numbers of unique protein sequences per barcode for patient 1
```
for f in *.aa.txt; do echo $f; awk '{if(NR==1){print};if($3=="productive"){split($2,n,":"); print n[2]":"n[3]":"$15,$0}}' $f | sort -u -k1,1 | wc -l; done
```
```
for f in *.aa.txt; do awk '$3=="productive"{split($2,n,":"); print n[2]":"n[3]":"$15,$0}' $f | sort -u -k1,1 | cut -f2- > ${f%.aa*}.unique.headless.txt; done
```