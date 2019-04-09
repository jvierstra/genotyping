#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-04
set ddir = ../data/$dstamp
set script = $HOME/Github/genotyping/pipeline/impute.sh
set vcfd = /net/seq/data/projects/altius-100-donors/archives/$dstamp/output
set vcfs = ($vcfd/filtered.all.hets-pass.recoded-final.vcf.gz $vcfd/filtered.all.hets-pass.recoded-final-ard.vcf.gz $vcfd/filtered.all.vcf.gz)
set baseoutd = ../results/$dstamp/dnaseI/output/phasing

source /net/module/Modules/default/tcsh

mkdir -p $ddir

foreach vcf ($vcfs)
  if ( ! -s $ddir/$vcf:t.tbi ) then
    module add htslib/1.7
    cp $vcf $ddir
    tabix -p vcf $ddir/$vcf:t
  endif

  set nm = $vcf:t:r:r
  set outd = $baseoutd/$nm

  bash $script $ddir/$vcf:t $outd
end

exit 0
