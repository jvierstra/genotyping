#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

# Richard says sequencing_depth for RNA and mapping_nuclear for DNaseI

set dna = production_qc_dnase_20190424-120705.parsed2.csv
set rna = production_qc_rna_20190424-120657.parsed2.csv

cut -f7,14 -d',' $rna \
  | tr ',' '\t' \
 >! rna.readcount

cut -f7,16 -d',' $dna \
  | tr ',' '\t' \
 >! dna.readcount

exit 0
