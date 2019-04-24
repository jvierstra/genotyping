#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-04
set ddir = `readlink -f ../data/$dstamp`
set baseind = ../results/$dstamp/rnaseq.filteredR2/output/imputed/filtered.all.renamed-cols
set outd = ../results.stats/rnaseq/$dstamp

mkdir -p $outd

set samplecounts = $ddir/prodQC/rna.readcount
set donorsmatched = $ddir/prodQC/donors.matched.simple
cp $samplecounts $outd
awk 'NR > 1' $samplecounts \
  | awk '{ s+=$2 } END { print int(s/NR) }' \
 >! $outd/samplecounts.avg.txt

awk 'BEGIN {OFS="\t"; print "Donor", "TissueCulture", "DNaseI-ID", "RNA-seq-ID", "Individual" } ; { print }' $donorsmatched \
 >! $outd/$donorsmatched:t

# haplotypes
mkdir -p $outd/by-individual
foreach t (haplotypes haplotypic_counts allelic_counts allele_config)
  foreach indiv (`find $baseind/ -maxdepth 1 -mindepth 1 -type d`)
    set fs = (`find $indiv/ -maxdepth 1 -mindepth 1 -type f -name "*.$t.txt" | sort`)
    awk 'NR == 1' $fs[1] \
     >! $outd/.header

    awk 'FNR > 1' $fs \
      | sort -k1,1 \
      | cat $outd/.header - \
     >! $outd/by-individual/$indiv:t.$t.txt

    rm -f $outd/.header
  end
end

# allelic expression


exit 0
