#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

# requires ../src.plots/ae-scatterplot.tcsh to have been run

set dstamp = 2019-04-24
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


# summarize FDR thresholded results
set baseind = ../results.plots/$dstamp/rnaseq-ae.filteredR2/output/expression/filtered.all

set fs = (`find $baseind/ -maxdepth 1 -mindepth 1 -type f -name "*.txt"`)
set fs = (`echo $fs | tr ' ' '\n' | awk '{ n=split($1, a, "/"); print $1"\t"substr(a[n], 7); }' | sort -nk2,2 | cut -f1`)

awk 'FNR > 1' $fs \
  | sort -sk4,4 --buffer-size=16G \
  | awk 'BEGIN {OFS="\t"} ; { \
          id = $1"|"$2"|"$3"|"$4; \
          counts[id]+=1; \
          trailing = substr($(NF-1),7); \
          if ( $7 > 0 ) { \
            pos[id]+=1; \
            if ( ("pos", id) in indiv ) { \
              trailing = indiv["pos", id]";"trailing; \
            } \
            indiv["pos", id]=trailing; \
          } else { \
            neg[id]+=1; \
            if ( ("neg", id) in indiv ) { \
              trailing = indiv["neg", id]";"trailing; \
            } \
            indiv["neg", id]=trailing; \
          } \
        } END { \
          for ( c in counts ) { \
            if ( !(c in pos) ) { \
              pos[c]=0; \
              indiv["pos", c] = "NA"; \
            } \
            if ( !(c in neg) ) { \
              neg[c]=0; \
              indiv["neg", c] = "NA"; \
            } \
            print c, pos[c], neg[c], pos[c]+neg[c], indiv["pos", c], indiv["neg", c]; \
          } \
        }' \
  | sort -srnk4,4 --buffer-size=16G \
  | tr '|' '\t' \
  | awk 'BEGIN { \
          OFS="\t"; \
          print "chr", "start", "stop", "gene-transcript", "+NumberIndiv", "-NumberIndiv", "TotalIndiv", "+IndivID", "-IndivID"; \
        } ; { print }' \
 >! $outd/ae-indiv-per-gene-fdr0.01.txt

exit 0
