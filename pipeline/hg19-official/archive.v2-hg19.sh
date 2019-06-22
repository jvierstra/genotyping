#!/bin/bash -x

#data_dir=$(dirname "$0")/../data
cd "$(dirname "$0")/.."

today=$(date -I)
archive=archives.no36bpfilter/$today.tar
dir=archives.no36bpfilter/$today

mkdir -p "$dir"

mkdir -p "$dir/input"
mkdir -p "$dir/output"

idir=data/input.no36bpfilter
odir=data/output.no36bpfilter

cp -r scripts2.hg19-no_36bp_filter/ "$dir/"
(cp $idir/* "$dir/input/") &
cp $odir/filtered.all.vcf.gz "$dir/output/"
cp $odir/filtered.all.vcf.gz.tbi "$dir/output/"
cp $odir/filtered.all.hets-pass.recoded-final.vcf.gz "$dir/output/"
cp $odir/filtered.all.hets-pass.recoded-final.vcf.gz.tbi "$dir/output/"
cp $odir/individual.counts/merged.all.vcf.gz "$dir/output/filtered.all.hets-pass.recoded-final-ard.vcf.gz"
cp $odir/individual.counts/merged.all.vcf.gz.tbi "$dir/output/filtered.all.hets-pass.recoded-final-ard.vcf.gz.tbi"
cp $odir/individual.counts/merged.all.ai.bed "$dir/output/filtered.all.hets-pass.recoded-final-ard.ai.bed"

module add samtools/1.7
module add htslib/1.7
module add vcftools

wait

tar -cf "$archive" "$dir"
