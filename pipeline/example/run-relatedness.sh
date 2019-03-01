#!/bin/bash -x

input=filtered.hwe.0.01.vcf.gz
output=relatedness.txt
makescript=../cluster_relatedness.sh

vcftools --gzvcf $input --relatedness --stdout > $output

out=final

# discard any -nan value (given when BAM is too sparse)
f=$output
any=$(awk '$NF == "-nan"' "$f" | awk 'END {print NR}')
if [[ $any -gt 0 ]] ; then
  awk '$NF != "-nan"' "$f" \
    > "$f.no-nans"
  f=$f.no-nans
fi

source "$makescript" "$f" "$out"

python3 cluster.py

exit 0
