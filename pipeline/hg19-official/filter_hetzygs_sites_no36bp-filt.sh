#!/bin/bash -x

usage () {
  echo -e "Usage: $0 <output-dir>" >&2
  exit 2
}

if [[ $# != 1 ]]; then
   usage
fi

output_dir=$1

module load bedops
module load bcftools/1.7
module load htslib/1.7
module load vcftools/0.1.14
module load python/2.7.11
#module load pysam


#gather all sites with at least 1 hetzygs indvidual and >30 reads total 
set -x
for i in ${output_dir}/samples/*.bed ; do 
  tmp=$output_dir/samples/sorted.$(basename "$i")
  sort-bed "$i" > "$tmp"
  mv "$tmp" "$i"
done
bedops --ec -u ${output_dir}/samples/AG*.bed | awk '$4!=$5' \
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
bedops -n 1 ${output_dir}/hetzygs_GQ50_DP.bed /net/seq/data/projects/altius-100-donors/hg19/data.references/hg19.blacklist.bed \
> ${output_dir}/hetzygs_GQ50_DP.no-ENCODE_blacklist.bed

#removed 36-bp filter

#Filter VCF file for SNVs that pass filters; note that this file is still intermediate
bcftools filter -R ${output_dir}/hetzygs_GQ50_DP.no-ENCODE_blacklist.bed -Oz ${output_dir}/filtered.all.vcf.gz > ${output_dir}/filtered.all.hets-pass.vcf.gz
tabix -p vcf ${output_dir}/filtered.all.hets-pass.vcf.gz

# bgzip and tabix up the individual genotype information
for sample in $(ls ${output_dir}/samples/*.bed); do bgzip -c $sample > ${sample}.gz; tabix -p bed ${sample}.gz; done

# recode filtered individual (sample) genotypes
python2  /net/seq/data/projects/altius-100-donors/scripts/genotyping/scripts/recode_vcf.py ${output_dir}/filtered.all.hets-pass.vcf.gz ${output_dir}/samples/ ${output_dir}/filtered.all.hets-pass.recoded.vcf.gz


# recompute VCF metrics such as allele frequency, etc.
vcftools --stdout --recode --recode-INFO-all --gzvcf ${output_dir}/filtered.all.hets-pass.recoded.vcf.gz | bgzip -c > ${output_dir}/filtered.all.hets-pass.recoded-final.vcf.gz
tabix -p vcf ${output_dir}/filtered.all.hets-pass.recoded-final.vcf.gz

#vcftools --out  ${output_dir}/filtered.all.hets-pass.recoded-final --gzvcf ${output_dir}/filtered.all.hets-pass.recoded-final.vcf.gz --relatedness --het --depth --TsTv-by-count --TsTv-by-qual

#zcat filtered.all.hets-pass.recoded-final.vcf.gz | perl -ne 'print "$1\n" if /AF1=([^,;]+)/'
