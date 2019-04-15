#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-04
set ddir = `readlink -f ../data/$dstamp`
set genotypeOutd = ../results.stats/dnaseI.genotypes/$dstamp
set imputeInd = ../results/$dstamp/dnaseI-SNV-combined/output/phasing-imputation
set imputeOutd = ../results.stats/dnaseI-SNV-combined/$dstamp

(workers/dnaseI.stats.imputation.tcsh $ddir $imputeInd $imputeOutd) &
(workers/dnaseI.stats.genotyping.tcsh $ddir $genotypeOutd) &

wait

exit 0
