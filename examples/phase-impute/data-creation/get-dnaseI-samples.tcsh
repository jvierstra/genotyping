#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

# It's a good idea to compare sample.order and sample.order.renamed against earlier dates
#  Right now, new samples seem to get appended due to increasing AG numbers

set dstamp = 2019-04-09
set vcf = /net/seq/data/projects/altius-100-donors/archives/$dstamp/output/filtered.all.vcf.gz
set outf = sample.order

source /net/module/Modules/default/tcsh
module add bcftools

bcftools query -l $vcf \
 >! sample.order

awk '{ print "Indiv-"NR }' sample.order \
 >! sample.order.renamed

exit 0
