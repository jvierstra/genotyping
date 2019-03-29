#!/bin/bash

output_dir=/net/seq/data/projects/genotyping/results.dgf-samples.merge2.genotype




#gather all sites with at least 1 hetzygs indvidual and >30 reads total 
bedops --ec -u ${output_dir}/samples/individual.*.bed | awk '$4!=$5' \
| awk -v OFS="\t" ' \
	function print_row() { print last_chrom, last_start, last_start+1, ".", dp, n_hets; } \
	{ \
		if (last_chrom!=$1 || last_start!=$2) { \
			if(NR>1) { print_row(); } \
			last_chrom=$1; last_start=$2; dp=$7; n_hets=1; \
		} else { \
			dp+=$7; n_hets+=1; \
		} \
	} \
	END { print_row(); }' \
| awk '$5>30' > ${output_dir}/hetzygs_GQ50_DP.bed

#remove positions that fall in ENCODE blacklist regions
cat ${output_dir}/hetzygs_GQ50_DP.bed \
| bedops -n -1 - /home/ehaugen/refseq/hg38/hg38.blacklist.bed \
> ${output_dir}/hetzygs_GQ50_DP.no-ENCODE_blacklist.bed

#Filter SNVs that are within 36 bp of eachother
closest-features --ec --delim '\t' --dist --closest --no-overlaps \
	${output_dir}/hetzygs_GQ50_DP.no-ENCODE_blacklist.bed ${output_dir}/hetzygs_GQ50_DP.no-ENCODE_blacklist.bed \
| awk '$(NF)>36 || $(NF)<-36' \
| cut -f1-6 > ${output_dir}/hetzygs_GQ50_DP.no-ENCODE_blacklist.36bp-window.bed

#Filter VCF file for SNVs that pass filters; note that this file is still intermediate
bcftools filter -R ${output_dir}/hetzygs_GQ50_DP.no-ENCODE_blacklist.36bp-window.bed -Oz ${output_dir}/filtered.all.vcf.gz > ${output_dir}/filtered.all.hets-pass.vcf.gz
tabix -p vcf ${output_dir}/filtered.all.hets-pass.vcf.gz

# bgzip and tabix up the individual genotype information
for sample in $(ls ${output_dir}/samples/*.bed); do bgzip -c $sample > ${sample}.gz; tabix -p bed ${sample}.gz; done

# recode filtered individual (sample) genotypes
python2 /home/jvierstra/proj/code/genotyping/scripts/recode_vcf.py ${output_dir}/filtered.all.hets-pass.vcf.gz ${output_dir}/samples/ ${output_dir}/filtered.all.hets-pass.recoded.vcf.gz

# recompute VCF metrics such as allele frequency, etc.
vcftools --stdout --recode --recode-INFO-all --gzvcf ${output_dir}/filtered.all.hets-pass.recoded.vcf.gz | bgzip -c > ${output_dir}/filtered.all.hets-pass.recoded-final.vcf.gz
tabix -p vcf ${output_dir}/filtered.all.hets-pass.recoded-final.vcf.gz

#vcftools --out  ${output_dir}/filtered.all.hets-pass.recoded-final --gzvcf ${output_dir}/filtered.all.hets-pass.recoded-final.vcf.gz --relatedness --het --depth --TsTv-by-count --TsTv-by-qual

#zcat filtered.all.hets-pass.recoded-final.vcf.gz | perl -ne 'print "$1\n" if /AF1=([^,;]+)/'
