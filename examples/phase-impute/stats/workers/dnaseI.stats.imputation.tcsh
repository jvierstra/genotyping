#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

if ( $#argv != 3 ) then
  printf "Need <data-dir> <input-dir> <base-outdir> args\n"
  exit -1
endif

set ddir = $1
set ind = $2
set outd = $3

set r2thresh = 0.3
set indirImpute = $ind
set genotypes = ($ddir/filtered.all.renamed-cols.vcf.gz)
set imputeOutd = $outd

source /net/module/Modules/default/tcsh
module add bcftools/1.7

foreach vcf ($genotypes)
  set vcfnm = $vcf:t:r:r

  set sampleorder = (`bcftools query -l $vcf`)

  ####################
  # Imputation
  ####################
  set indir = $indirImpute/$vcfnm
  set outdir = $imputeOutd/$vcfnm
  mkdir -p $outdir

  # number of imputed positions with R2 >= $r2thresh or ground truth SNVs that imputation did not include
  foreach chr (`find $indir/ -maxdepth 1 -mindepth 1 -type f -name "*.vcf.gz"`)
    set chrnm = `echo $chr:t | cut -f1 -d'.'`
    (bcftools view -i 'R2>='$r2thresh' || VDB!="."' $chr \
      | awk '$1 !~ /^#/' \
      | awk -v c=$chrnm 'END { print c"\t"NR; }' \
     >! $outdir/$vcfnm.imputed.R2.$chrnm.details) &
  end

  # number of individuals who are het per imputed site
  foreach chr (`find $indir/ -maxdepth 1 -mindepth 1 -type f -name "*.vcf.gz"`)
    set chrnm = `echo $chr:t | cut -f1 -d'.'`
    (bcftools view $chr \
      | awk -F"\t" '$1 !~ /^#/ {print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t" gsub(/0\|1|1\|0|0\/1|1\/0/,"")}' \
      | awk '$NF > 0' \
     >! $outdir/$vcfnm.$chrnm.heterozygous.nonzero.perImputed) &
  end

  # number of hets per individual
  foreach chr (`find $indir/ -maxdepth 1 -mindepth 1 -type f -name "*.vcf.gz"`)
    set chrnm = `echo $chr:t | cut -f1 -d'.'`
    foreach indiv ($sampleorder)
      (bcftools view -s $indiv $chr \
        | awk -F"\t" '$1 !~ /^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t" gsub(/0\|1|1\|0|0\/1|1\/0/,"")}' \
        | awk '$NF > 0' \
        | awk -v i=$indiv -v c=$chrnm 'END { print i"\t"c"\t"NR }' \
       >! $outdir/$vcfnm.$indiv.$chrnm.nhets) &
    end
    wait
  end

  wait

  rm -f $outdir/$vcfnm.indiv.nhets.summary
  foreach indiv ($sampleorder)
    set fs = (`find $outdir/ -maxdepth 1 -mindepth 1 -type f -name "$vcfnm.$indiv.*.nhets"`)
    cat $fs \
      | awk -v i=$indiv '{ s+=$NF } END { print i"\t"s }' \
    >>! $outdir/$vcfnm.indiv.nhets.summary

    foreach f ($fs)
      rm -f $f
    end
  end

  set fs = (`find $outdir/ -maxdepth 1 -mindepth 1 -type f -name "$vcfnm.imputed.R2.*.details"`)
  cat $fs \
    | awk '{ s+=$2 } END { print s }' \
   >! $outdir/$vcfnm.imputed.R2.summary

  foreach f ($fs)
    rm -f $f
  end

  set fs = (`find $outdir/ -maxdepth 1 -mindepth 1 -type f -name "$vcfnm.*.heterozygous.nonzero.perImputed"`)
  sort-bed $fs \
    | cut -f1-2,4- \
   >! $outdir/$vcfnm.nindiv.heterozygous.nonzero.perImputed

  foreach f ($fs)
    rm -f $f
  end
end

exit 0
