#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-04
set r2thold = 0.3
set vcfd = ../data/$dstamp
set vcfs = ($vcfd/filtered.all.hets-pass.recoded-final.vcf.gz $vcfd/filtered.all.hets-pass.recoded-final-ard.vcf.gz $vcfd/filtered.all.vcf.gz)
set baseind = ../results/$dstamp/dnaseI/output/phasing-imputation
set baseoutd = ../results/$dstamp/dnaseI-SNV-combined/output/phasing-imputation

source /net/module/Modules/default/tcsh
module add bcftools/1.7
module add htslib/1.7
module add vcftools

foreach vcf ($vcfs)
  set vcfnm = $vcf:t:r:r
  set outd = $baseoutd/$vcfnm
  mkdir -p $outd

  foreach chr (`find $baseind/$vcfnm/ -maxdepth 1 -mindepth 1 -type f -name "*.vcf.gz"`)
    set chrnm = `echo $chr:t | cut -f1 -d'.'`

    # isec takes '-' but a separate tabix index is required...
    (bcftools view -i "R2 >= $r2thold" -Oz $chr \
      >! $outd/$chrnm.tmp1.vcf.gz && \
     tabix -p vcf $outd/$chrnm.tmp1.vcf.gz && \
     bcftools isec --complement $vcf $outd/$chrnm.tmp1.vcf.gz --regions $chrnm -w1 -Oz \
      >! $outd/$chrnm.tmp.vcf.gz && \
     tabix -p vcf $outd/$chrnm.tmp.vcf.gz && \
     bcftools concat $outd/$chrnm.tmp.vcf.gz $outd/$chrnm.tmp1.vcf.gz \
       | vcf-sort -p 4 \
       | bgzip -c \
      >! $outd/$chrnm.imputed.vcf.gz && \
     tabix -p vcf $outd/$chrnm.imputed.vcf.gz && \
     rm -f $outd/$chrnm.tmp.vcf.gz $outd/$chrnm.tmp1.vcf.gz && \
     rm -f $outd/$chrnm.tmp.vcf.gz.tbi $outd/$chrnm.tmp1.vcf.gz.tbi) &
  end
end

exit 0
