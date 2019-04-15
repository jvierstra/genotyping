#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

# Richard says sequencing_depth for RNA and mapping_nuclear for DNaseI

set dna = production_qc_dnase_20190405-172400.csv
set rna = production_qc_rna_20190405-172238.csv

cut -f9,18 -d',' $rna \
  | tr ',' '\t' \
 >! rna.readcount

cut -f9,20 -d',' $dna \
  | tr ',' '\t' \
 >! dna.readcount

exit 0
