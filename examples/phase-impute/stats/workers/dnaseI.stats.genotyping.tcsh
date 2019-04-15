#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

if ( $#argv != 2 ) then
  printf "Need <data-dir> <base-outdir> args\n"
  exit -1
endif

set ddir = $1
set outd = $2

set indirGenotype = $ddir
set genotypes = ($indirGenotype/filtered.all.vcf.gz)
set samplecounts = $ddir/prodQC/dna.readcount
set sampleorder_newnames = $ddir/sample.order.renamed
set sampleorder_oldnames = $ddir/sample.order

source /net/module/Modules/default/tcsh
module add bcftools/1.7

mkdir -p $outd
paste $sampleorder_newnames $sampleorder_oldnames \
 >! $outd/sample.ordering.txt
cp $samplecounts $outd
awk 'NR > 1' $samplecounts \
  | awk '{ s+=$2 } END { print int(s/NR) }' \
 >! $outd/samplecounts.avg.txt

# genotype stats
foreach vcf ($genotypes)
  set vcfnm = $vcf:t:r:r
  set outdir = $outd
  mkdir -p $outdir

  ####################
  # Genotype
  ####################
  # how many SNVs per individual?
  (zcat $vcf \
    | awk '$1 !~ /^##/' \
    | cut -f10- \
    | awk 'BEGIN {OFS="\t"} ; { \
            if ( NR > 1 ) { \
              for (i=1; i<=NF; ++i) { \
                split($i, a, ":"); \
                v=0; \
                if ( a[3] > 0 ) { v=1; } \
                rows[i] += v; \
              } \
            } else { \
              mnf = NF; \
              print; \
            } \
          } END { \
            printf "%s", rows[1]; \
            for(i=2; i<=mnf; ++i) { \
              printf "\t%s", rows[i]; \
            } \
            printf "\n"; \
          }' \
       >! $outdir/$vcfnm.snvs.counts && \
   awk 'NR == 2' $outdir/$vcfnm.snvs.counts \
     | awk '{ for(i=1;i<=NF;++i) { s+=$i; } t=NF } END { print int(s/t) }' \
    >! $outdir/$vcfnm.snvs.counts-average) &

  # number of individuals who are heterozygous per SNV
  (bcftools view $vcf \
    | awk -F"\t" '$1 !~ /^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t" gsub(/0\|1|1\|0|0\/1|1\/0/,"")}' \
    | tee $outdir/$vcfnm.nindiv.heterozygous.perSNV \
    | awk '$NF > 0' \
   >! $outdir/$vcfnm.nindiv.heterozygous.nonzero.perSNV) &

  # number of SNVs
  (zcat $vcf \
    | awk '$1 !~ /^#/' \
    | awk 'END { print NR }' \
   >! $outdir/$vcfnm.number-SNVs.txt) &

  # number of hets per individual
  set sampleorder = (`bcftools query -l $vcf`)
  foreach indiv ($sampleorder)
    (bcftools view -s $indiv $vcf \
      | awk -F"\t" '$1 !~ /^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t" gsub(/0\|1|1\|0|0\/1|1\/0/,"")}' \
      | awk '$NF > 0' \
      | awk -v i=$indiv 'END { print i"\t"NR }' \
     >! $outdir/$vcfnm.$indiv.nhets) &
  end

  wait

  rm -f $outdir/$vcfnm.indiv.nhets.summary
  foreach indiv ($sampleorder)
    set f = `find $outdir/ -maxdepth 1 -mindepth 1 -type f -name "$vcfnm.$indiv.nhets"`
    cat $f \
    >>! $outdir/$vcfnm.indiv.nhets.summary

    rm -f $f
  end
end

wait

exit 0
