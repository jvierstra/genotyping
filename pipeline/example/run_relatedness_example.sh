#!/bin/bash -x

input=filtered.hwe.0.01.vcf.gz
relative_thold=0.3
outdir=results.relatedness
script=../run_relatedness.sh

$script $input $relative_thold $outdir

exit 0
