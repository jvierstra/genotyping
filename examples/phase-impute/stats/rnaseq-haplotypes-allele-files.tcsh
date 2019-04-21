#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-04
set ddir = `readlink -f ../data/$dstamp`
set baseind = ../results.stats/rnaseq/$dstamp/by-individual
set outd = ../results.stats/rnaseq/$dstamp

set tmpdir = /tmp/`whoami`/$$

mkdir -p $outd
rm -rf $tmpdir
mkdir -p $tmpdir

set haplos = (`find $baseind/ -maxdepth 1 -mindepth 1 -type f -name "*.haplotypic_counts.txt"`)
set alleles = (`find $baseind/ -maxdepth 1 -mindepth 1 -type f -name "*.allelic_counts.txt"`)

awk 'BEGIN {OFS="\t"} ; { n=split(FILENAME, a, "/"); split(a[n], b, "."); print $1, $2, $2+1, $4"|"$5"|"b[1]"|"$6"|"$7; }' $alleles \
  | awk '$1 != "contig"' \
  | sort-bed - \
 >! $tmpdir/alleles

cut -f1-3 $tmpdir/alleles \
  | uniq \
  | bedmap --exact --delim "\t" --echo --echo-map-id --count - $tmpdir/alleles \
 >! $outd/allele_counts.bed

set hs = ()
foreach h ($haplos)
  tr '\t' '?' < $h \
   >! $tmpdir/$h:t

  set hs = ($hs $tmpdir/$h:t)
end


awk 'BEGIN {FS="?"; OFS="\t"} ; { n=split(FILENAME, a, "/"); split(a[n], b, "."); print $1, $2, $2+1, $8"|"$9"|"b[1]"|"$10"|"$11; }' $hs \
  | awk '$1 != "contig"' \
  | sort-bed - \
 >! $tmpdir/haplos

cut -f1-3 $tmpdir/haplos \
  | uniq \
  | bedmap --exact --delim "\t" --echo --echo-map-id --count - $tmpdir/haplos \
 >! $outd/haplotypes_counts.bed

rm -rf $tmpdir

exit 0
