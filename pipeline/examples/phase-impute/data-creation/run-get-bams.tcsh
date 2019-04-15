#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set tagid = genotyping-100-donors-t-cells-rna-datasets

./get-bams.bash $tagid bamlist.txt

mkdir -p bams

foreach b (`cat bamlist.txt`)
  set nm = `echo $b | cut -f6,7 -d'/' | tr '/' '-' | sed 's;aggregation-;AG;'`
  cp -s $b bams/$nm.bam
  cp -s $b.bai bams/$nm.bam.bai
end

exit 0
