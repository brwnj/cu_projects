#poly(a)

getting the counts
```
bedtools map -c 4 -o max -null 0 -a polya_db_slop.bed.gz -b MP55.pos.bedgraph.gz | cut -f 4,7
```
generating "replicates"
```
for f in *.counts;\
    do awk 'BEGIN{OFS=FS="\t"}{foo=int(rand()*10);if($2!=0){print $1,$2+foo}else{print}}' $f\
    > ${f%.counts}.x.counts;\
done
```