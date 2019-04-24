#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-04
set ddir = `readlink -f ../data/$dstamp`
set baseind = ../results/$dstamp/rnaseq-ae.filteredR2/output/expression/filtered.all.renamed-cols
set input = $baseind/final.ae.best.with-header.uniqs.filtered.txt
set base_tlabel = "log2(foldchange)"
set threshold_label = "abs("$base_tlabel")>=1"
set outd = ../results.stats/rnaseq.ae/$dstamp

mkdir -p $outd

# number of genes per individual
awk 'BEGIN {FS="\t";OFS="\t"} ; { \
      n=split($4,a,";"); \
      print a[1]";"a[2], $7, $9; \
    }' $input \
  | awk -v l=$threshold_label \
      'BEGIN {OFS="\t"} ; { \
        if ( NR > 1 ) { \
          indiv[$NF]+=1; \
        } \
      } END { \
        for ( i in indiv ) { \
          print i, indiv[i]; \
        } \
      }' \
  | awk '{ $1=substr($1,7); print }' \
  | sort -nk1,1 \
  | awk '{ print "Indiv-"$1, $2 }' \
  | awk -v l=$threshold_label \
        'BEGIN { \
          OFS="\t"; \
          label="#allelic-expressed_"l; \
          print "Indiv", label; \
        } ; { \
          print; \
        }' \
 >! $outd/ae-genes-per-indiv.txt

# indviduals per gene
awk 'BEGIN {OFS="\t"} ; { \
      if ( NR > 1 ) { \
        if ( $7 >= 0 ) { \
          genesplus[$4]+=1; \
        } else { \
          genesminus[$4]+=1; \
        } \
      } \
    } END { \
      for ( g in genesplus ) { \
        if ( ! (g in genesminus) ) { \
          genesminus[g]=0; \
        } \
        print g, genesplus[g], genesminus[g]; \
      } \
      for ( g in genesminus ) { \
        if ( ! (g in genesplus) ) { \
          genesplus[g]=0; \
          print g, genesplus[g], genesminus[g]; \
        } \
      } \
    }' $input \
  | sort -k1,1 \
  | awk -v l=$base_tlabel \
        'BEGIN { \
          OFS="\t"; \
          label="#Indiv_"l""; \
          print "Gene", "+_"label">=1", "-_"label"<=-1"; \
        } ; { \
          print; \
        }' \
 >! $outd/ae-indiv-per-gene.txt

exit 0
