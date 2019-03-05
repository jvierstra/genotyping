#!/bin/bash -x

cluster_file=filtered.hwe.0.01.vcf.relatedness.greedy # produced from ../run_relatedness.sh
output_dir=results.merged
script=../merge_bams_by_individual.sh
mkdir -p $output_dir

chmod 755 $script
$script $cluster_file $output_dir

exit 0
