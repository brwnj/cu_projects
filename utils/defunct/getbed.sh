mysql5 --user=genome --host=genome-mysql.cse.ucsc.edu -A -D mm9 -P 3306 -e "select chrom, txStart, txEnd, E.name, X.value as geneName, strand, exonStarts, exonEnds from ensGene as E, ensemblToGeneName as X where X.name=E.name;" > tmp.notbed 
 awk 'BEGIN { OFS = "\t"; } ;
    (NR > 1){
            split($1, chr, /chr/)
            split($7, astarts, /,/);
            split($8, aends, /,/);
            starts=""
            sizes=""
            exonCount=0
            for(i in astarts){
                if (! astarts[i]) continue
                sizes=sizes""(aends[i] - astarts[i])","
                starts=starts""(astarts[i] = astarts[i] - $2)","
                exonCount=exonCount + 1
            }
            print chr[2],$2,$3,$4","$5,1,$6,$2,$3,".",exonCount,sizes,starts}'
            tmp.notbed 
 | sort -k1,1 -k2,2n > mm9.ens.bed
 rm tmp.notbed
 
 
 mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D $ORG -P 3306   -e "select chrom,txStart,txEnd,K.name,X.geneSymbol,strand,exonStarts,exonEnds from knownGene as K,kgXref as X where  X.kgId=K.name;" > tmp.notbed
 grep -v txStart tmp.notbed | awk '
         BEGIN { OFS = "\t"; FS = "\t"} ;
             {
                 delete astarts;
                 delete aends;
                 split($7, astarts, /,/);
                 split($8, aends, /,/);
                 starts=""
                 sizes=""
                 exonCount=0
                 for(i=1; i <= length(astarts); i++){
                     if (! astarts[i]) continue
                     sizes=sizes""(aends[i] - astarts[i])","
                     starts=starts""(astarts[i] = astarts[i] - $2)","
                     exonCount=exonCount + 1
                 }
                 print $1,$2,$3,$5","$4,1,$6,$2,$3,"0",exonCount,sizes,starts
             }' | sort -k1,1 -k2,2n > knownGene.$ORG.bed12