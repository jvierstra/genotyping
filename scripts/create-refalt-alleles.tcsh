#!/bin/tcsh -efx

# use topmed's ref/alt calls over the (large majority of) and map on individual counts
#  this will be the matrix used in count_genotypes.py
# an alternative could be to do a majority-wins ref/alt over our individuals
# $ind_pre_mtx is the output from pipeline/counts_tags_by_individual.sh

set topmedref_d = ../results.analysis.merged/results.cross.topmed/bed
set topmedrefs = ($topmedref_d/filtered.hwe.0.01.vcf.topmed.ALL.pass.bed $topmedref_d/filtered.hwe.0.01.vcf.topmed.ALL.TOPMed_freeze5_hg38_dbSNP.bed)
set ind_pre_mtx = ../results.analysis.merged/results.individual-counts/merged.counts.no-filter.bed
set outdir = ../results.analysis.merged/results.genotype.inmtx
mkdir -p $outdir

foreach t ($topmedrefs)
  set nm = $t:t.mtx

  (awk -F"\t" 'BEGIN {OFS="\t"} ; { print $1, $2, $3, $6"/"$7; }' $t \
    | bedmap --skip-unmapped --exact --echo --echo-map --delim "\t" - $ind_pre_mtx \
    | cut -f1-4,8- \
   >! $outdir/$nm) &
end

exit 0
