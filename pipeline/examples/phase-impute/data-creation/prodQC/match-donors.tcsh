#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dnase = production_qc_dnase_20190405-172400.csv
set rnase = production_qc_rna_20190405-172238.csv
set sampleorder = ../sample.order
set newsamplenames = ../sample.order.renamed

set donorcol = 16

setenv LANG C

awk 'NR > 1' $dnase \
  | tr ' ' '_' \
  | awk -F"," 'BEGIN {OFS="\t"} ; { for(i=1;i<=NF;++i) { if ($i == "") { $i= "NA" } } print }' \
  | sort -k16,16 \
 >! dnase.tmp

awk 'NR > 1' $rnase \
  | tr ' ' '_' \
  | awk -F"," 'BEGIN {OFS="\t"} ; { for(i=1;i<=NF;++i) { if ($i == "") { $i= "NA" } } print }' \
  | sort -k16,16 \
 >! rnase.tmp

join --check-order -1 16 -2 16 dnase.tmp rnase.tmp \
  | tr ' ' '\t' \
 >! donors.matched

cut -f1-2,10,47 donors.matched \
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
rm -f rnase.tmp

exit 0
