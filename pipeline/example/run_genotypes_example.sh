#!/bin/bash -x

###
# Example for hg38 TCells

ddir=../data/links # to all BAM files with corresponding index files
outdir=results.tcells
FASTA_CHROM_FILE=/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.chrom_sizes.bed
FASTA_FILE=/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.fa

script=../../call_genotypes.sh

chmod 755 $script
$script $ddir $outdir $FASTA_CHROM_FILE $FASTA_FILE

exit 0
