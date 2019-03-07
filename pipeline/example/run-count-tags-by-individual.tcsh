#!/bin/tcsh -efx
# author : sjn
# date : Mar.2019

set iter = iter1
set FASTA_CHROM_FILE = /net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.chrom_sizes.bed
set FASTA_FILE = /net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.fa
set GZVCF_FILE = ../results.tcells/$iter/filtered.hwe.0.01.vcf.gz

set output_dir = ../results.individual-counts/$iter
mkdir -p $output_dir/logs
cp ../../scripts/count_tags.py $output_dir

./counts_tags_by_individual.sh $FASTA_CHROM_FILE $FASTA_FILE $GZVCF_FILE $output_dir

exit 0
