#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dnase = production_qc_dnase_20190424-120705.parsed2.csv
set rnaseq = production_qc_rna_20190424-120657.parsed2.csv
set sampleorder = ../sample.order
set newsamplenames = ../sample.order.renamed

# donor coloums
set ddonorcol = 12
set rdonorcol = 12

# agg ID columns after joining things
set dag = 8
set rag = 39

setenv LANG C

awk 'NR > 1' $dnase \
  | tr ' ' '_' \
  | awk -F"," 'BEGIN {OFS="\t"} ; { for(i=1;i<=NF;++i) { if ($i == "") { $i= "NA" } } print }' \
  | tr ',' '\t' \
  | sort -k$ddonorcol,$ddonorcol \
 >! dnase.tmp

awk 'NR > 1' $rnaseq \
  | tr ' ' '_' \
  | awk -F"," 'BEGIN {OFS="\t"} ; { for(i=1;i<=NF;++i) { if ($i == "") { $i= "NA" } } print }' \
  | tr ',' '\t' \
  | sort -k$rdonorcol,$rdonorcol \
 >! rnaseq.tmp

join -i --check-order -1 $ddonorcol -2 $rdonorcol dnase.tmp rnaseq.tmp \
  | tr ' ' '\t' \
 >! donors.matched

cut -f1-2,$dag,$rag donors.matched \
 >! donors.matched.simple

sort -k3,3 donors.matched.simple \
 >! .tmp1

sed 's;\(.*\)\(/.*\).bam;\2;;' < $sampleorder \
  | sed 's;/;;;' \
  | paste - $newsamplenames \
  | join --check-order -1 3 -2 1 .tmp1 - \
  | awk '{ print $2"\t"$3"\t"$1"\t"$4"\t"$5; }' \
 >! .tmp2

mv .tmp2 donors.matched.simple

rm -f .tmp1
rm -f dnase.tmp
rm -f rnaseq.tmp

exit 0
