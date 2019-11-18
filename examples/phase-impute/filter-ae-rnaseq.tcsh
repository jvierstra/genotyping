#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-24
set baseind = ../results/$dstamp/rnaseq-ae.filteredR2/output/expression/filtered.all.renamed-cols
set calls = $baseind/final.ae.best.with-header.bed # already filtered out "inf"
set threshold = 1
set rosetta = ../data/$dstamp/prodQC/donors.matched.simple

awk 'NR > 1' $calls \
  | awk 'BEGIN { \
      OFS="\t" \
    } ; { \
      for ( i in c ) { delete c[i]; } \
      n=split($4, a, ";"); \
      for(i=1;i<=n;++i) { \
        split(a[i], b, ":"); \
        c[b[1]] = b[2]; \
      } \
      id=c["gene_id"]";"c["gene_name"]";"c["gene_type"]; \
      if ( "transcript_id" in c ) { \
        id=id";"c["transcript_id"]; \
      } \
      print $1, $2, $3, id, $5, $6, $8, $9, $12; \
    }' \
  | uniq \
  | sort --buffer-size=8G \
  | uniq \
  | awk -v r=$rosetta 'BEGIN { \
      OFS="\t"; \
      print "chr", "start", "stop", "gene", "refCount", "altCount", "log2_aFC", "n_variants", "sample-id"; \
      while ( (getline line < r) > 0 ) { \
        split(line, c, "\t"); \
        names[c[4]]=c[5]; \
      } \
    } ; { \
      sanity=10; \
      minr=5; \
      if ( $5>=$6*sanity || $6>=$5*sanity ) { \
        $1="***"$1; \
      } else if ( $5<=minr || $6<=minr ) { \
        $1="***"$1; \
      } \
      split($NF, a, "-"); \
      ag=a[2]; \
      $NF = names[ag]; \
      print; \
    }' \
  | tee $baseind/$calls:t:r.uniqs.txt \
  | awk '$1 !~ /^\*/' \
  | sort -nk7,7 \
  | (head -500; tail -500) \
 >! $baseind/$calls:t:r:r.tophits.filtered.txt.tmp

awk '$1 !~ /^*/' $baseind/$calls:t:r.uniqs.txt \
  | awk -v t=$threshold '$7 >= t || $7 <= -t' \
 >! $baseind/$calls:t:r.uniqs.filtered.txt

awk 'NR == 1' $baseind/$calls:t:r.uniqs.txt \
  | cat - $baseind/$calls:t:r:r.tophits.filtered.txt.tmp \
 >! $baseind/$calls:t:r:r.tophits.filtered.txt

rm -f $baseind/$calls:t:r:r.tophits.filtered.txt.tmp

exit 0
