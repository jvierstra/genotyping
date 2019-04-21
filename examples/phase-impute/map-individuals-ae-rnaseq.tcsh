#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-04
set baseind = ../results/$dstamp/rnaseq-ae.filteredR2/output/expression/filtered.all.renamed-cols
set calls = $baseind/final.ae.best.with-header.bed # already filtered out "inf"
set rosetta = ../data/$dstamp/prodQC/donors.matched.simple
set baseoutd = ../results.simple/$dstamp/allelic.expression.renamed-cols

mkdir -p $baseoutd

if ( ! -s $baseoutd/final.ae.premap.bed ) then
awk 'NR > 1' $calls \
  | awk 'BEGIN { \
      OFS="\t" \
    } ; { \
      n=split($4, a, ";"); \
      for(i=1;i<=n;++i) { \
        split(a[i], b, ":"); \
        c[b[1]] = b[2]; \
      } \
      nm = c["gene_id"]";"c["gene_name"]; \
      if ( "transcript_id" in c ) { \
        nm = nm";"c["transcript_id"]; \
      } \
      print $1, $2, $3, nm";"$12, $8, $5, $6; \
    }' \
  | awk -v r=$rosetta 'BEGIN { \
      OFS="\t"; \
      while ( (getline line < r) > 0 ) { \
        split(line, c, "\t"); \
        names[c[4]]=c[5]; \
      } \
    } ; { \
      n=split($4, a, ";"); \
      lnag=a[n]; \
      split(lnag, ag, "-"); \
      id = names[ag[2]]; \
      $4 = a[1]; \
      for(i=2;i<n;++i) { \
        $4 = $4";"a[i]; \
      } \
      sanity=10; \
      minr=5; \
      marker=""; \
      if ( $(NF-1)>=$NF*sanity || $NF>=$(NF-1)*sanity ) { \
        marker="***"; \
      } else if ( $(NF-1)<=minr || $NF<=minr ) { \
        marker="^^^"; \
      } \
      $4=$4";"id"?"$5"?"$(NF-1)"?"$NF; \
      print $0; \
    }' \
  | sort-bed - \
 >! $baseoutd/final.ae.premap.bed
endif

cut -f1-4 $baseoutd/final.ae.premap.bed \
  | awk 'BEGIN {OFS="\t"} ; { \
          n=split($4, a, ";"); \
          idv=a[n]; \
          split(idv, details, "?"); \
          idv=details[1]; \
          score=details[2]; \
          refcnt=details[3]; \
          altcnt=details[4]; \
          geneinfo=""; \
          for(i=1;i<n;++i) { \
            geneinfo=geneinfo";"a[i]; \
          } \
          print $1, $2, $3, idv":"score, substr(geneinfo, 2); \
        }' \
  | bedmap --exact --echo --echo-map-id-uniq --delim "\t" - \
  | cut -f1-3,5- \
  | uniq \
  | sed 's;Indiv-;;g' \
 >! $baseoutd/final.ae.bed

exit 0
