#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-24
set ddir = ../data/$dstamp
set script = ~sjn/Github/genotyping/pipeline/impute.sh
set vcfd = /net/seq/data/projects/altius-100-donors/archives/$dstamp/output
#set vcfs = ($vcfd/filtered.all.vcf.gz $vcfd/filtered.all.hets-pass.recoded-final.vcf.gz $vcfd/filtered.all.hets-pass.recoded-final-ard.vcf.gz)
set vcfs = ($vcfd/filtered.all.vcf.gz)
set baseoutd = ../results/$dstamp/dnaseI/output/phasing-imputation

source /net/module/Modules/default/tcsh
module add bcftools/1.7
setenv LANG C

mkdir -p $ddir

foreach vcf ($vcfs)
  set newsamplenames = `readlink -f $ddir/sample.order.renamed`

  set vcf_renamed = $ddir/$vcf:t:r:r.renamed-cols.vcf.gz

  if ( ! -s $vcf_renamed ) then
    module add htslib/1.7

    bcftools reheader --samples $newsamplenames $vcf \
     >! $vcf_renamed

    tabix -p vcf $vcf_renamed
  endif

  set nm = $vcf_renamed:t:r:r
  set outd = $baseoutd/$vcf_renamed:t:r:r
  mkdir -p $outd

  bash $script $vcf_renamed $outd
end

exit 0
