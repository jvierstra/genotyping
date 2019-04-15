#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-04
set ddir = ../data/$dstamp
set script = ~sjn/Github/genotyping/pipeline/impute.sh
set vcfd = /net/seq/data/projects/altius-100-donors/archives/$dstamp/output
set vcfs = ($vcfd/filtered.all.hets-pass.recoded-final.vcf.gz $vcfd/filtered.all.hets-pass.recoded-final-ard.vcf.gz $vcfd/filtered.all.vcf.gz)
set baseoutd = ../results/$dstamp/dnaseI/output/phasing-imputation

source /net/module/Modules/default/tcsh
module add bcftools/1.7
setenv LANG C

mkdir -p $ddir

foreach vcf ($vcfs)
  set sampleorder = `readlink -f $ddir/sample.order`
  set newsamplenames = `readlink -f $ddir/sample.order.renamed`

  set nm = $vcf:t:r:r
  set outd = $baseoutd/$nm
  mkdir -p $outd

  bcftools query -l $vcf \
    | sort \
   >! $sampleorder

  awk '{ print "Indiv-"NR }' $sampleorder \
   >! $newsamplenames

  if ( ! -s $ddir/$vcf:t.tbi ) then
    module add htslib/1.7
    bcftools annotate --samples-file $sampleorder $vcf \
     | bcftools reheader --samples $newsamplenames \
     | bgzip -c > $ddir/$vcf:t

    tabix -p vcf $ddir/$vcf:t
  endif

  bash $script $ddir/$vcf:t $sampleorder $newsamplenames $outd
end

exit 0
