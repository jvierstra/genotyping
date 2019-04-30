#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-24
set ddir = `readlink -f ../data/$dstamp`
set pipelineAI = /net/seq/data/projects/altius-100-donors/archives/$dstamp/output/filtered.all.hets-pass.recoded-final-ard.ai.bed
set genotypeOutd = ../results.stats/dnaseI.genotypes/$dstamp
set imputeInd = ../results/$dstamp/dnaseI-SNV-combined/output/phasing-imputation
set imputeOutd = ../results.stats/dnaseI-SNV-combined/$dstamp

mkdir -p $genotypeOutd
awk '$NF != "nan"' $pipelineAI \
  | sort-bed - \
  | tee $genotypeOutd/$pipelineAI:t \
  | awk 'END { print NR }' \
 >! $genotypeOutd/$pipelineAI:t:r.counts

(workers/dnaseI.stats.imputation.tcsh $ddir $imputeInd $imputeOutd) &
(workers/dnaseI.stats.genotyping.tcsh $ddir $genotypeOutd) &

wait

exit 0
