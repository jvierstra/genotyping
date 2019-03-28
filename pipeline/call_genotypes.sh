#!/bin/bash

###
# TODO:
#	-Check relatedness (vcf-tools)
#	-Filter SNPs
#

FASTA_CHROM_FILE=/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.chrom_sizes.bed
FASTA_FILE=/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.fa

output_dir=/net/seq/data/projects/genotyping/results.dgf-samples.merge2.genotype

# Make BAM file list
#cut -f3  /home/jvierstra/proj/ftd.updated/hg38.DGF.ver4.params.txt > ${output_dir}/filelist.txt
ls /net/seq/data/projects/genotyping/results.dgf-samples.merge2/*.bam > ${output_dir}/filelist.txt

#params:
min_Q=500
min_GQ=50
min_DP=12
min_AD=4
hwe_cutoff=0.01

###DO NOT EDIT BELOW
rm -rf ${output_dir}/logs && mkdir -p ${output_dir}/logs ${output_dir}/samples

# Make chunks
cat ${FASTA_CHROM_FILE} \
| grep -v random | grep -v _hap | grep -v chrM | grep -v chrUn | grep -v chrX | grep -v chrY | grep -v chrEBV \
| cut -f1,3 \
| awk -v step=5000000 -v OFS="\t" '{ for(i=step; i<=$2; i+=step) { print $1":"i-step+1"-"i; } print $1":"i-step+1"-"$2; }' \
> ${output_dir}/inputs.txt

njobs=$(wc -l < ${output_dir}/inputs.txt)


cat <<__SCRIPT__ > ${output_dir}/slurm.bam_call_genotypes_chunk
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%A.%a.out
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=queue0

set -e -o pipefail

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

module load samtools/1.7
module load bcftools/1.7
module load vcftools/0.1.14
module load htslib/1.7

region=\`cat ${output_dir}/inputs.txt | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1\`

bcftools mpileup -r \${region} \
	-Q 20 -d 1000 -I -E \
	-f ${FASTA_FILE} \
	-b ${output_dir}/filelist.txt \
	-a FORMAT/DP,FORMAT/AD \
| bcftools call -f GQ -cv -Ov \
| vcftools --stdout \
	--minQ ${min_Q} --minGQ ${min_GQ} --minDP ${min_DP} --max-alleles 2 --hwe ${hwe_cutoff}\
	--recode --recode-INFO-all --vcf - \
| bgzip -c > ${output_dir}/\${region}.filtered.vcf.gz

tabix -p vcf ${output_dir}/\${region}.filtered.vcf.gz

rm -rf \${TMPDIR}
__SCRIPT__

cat <<__SCRIPT__ > ${output_dir}/slurm.bam_call_genotypes_merge
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%J.out
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=queue0

set -euxo pipefail

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

module load bcftools/1.7
module load vcftools/0.1.14
module load htslib/1.7

cat inputs.txt | awk -v dir="${output_dir}" '{ print dir"/"\$1".filtered.vcf.gz"; }' > \${TMPDIR}/mergelist.txt

bcftools concat -f \${TMPDIR}/mergelist.txt -Ou \
| bcftools annotate -Oz -a /home/jvierstra/data/dbSNP/v151/All_20180418.fixed-chrom.vcf.gz -c ID --threads 16 \
> ${output_dir}/filtered.all.vcf.gz
tabix -p vcf ${output_dir}/filtered.all.vcf.gz

python2 ${script_dir}/filter_by_sample.py \
	--min_dp ${min_DP}  --min_ad ${min_AD} --min_gq ${min_GQ} \
	${output_dir}/filtered.all.vcf.gz ${output_dir}/samples \
> ${output_dir}/filtered.all.bed

__SCRIPT__

JOB0=$(sbatch --export=ALL \
	--job-name=genotyping \
	--array=1-${njobs} \
	${output_dir}/slurm.bam_call_genotypes_chunk)
echo $JOB0

JOB1=$(sbatch --export=ALL \
	--job-name=genotyping.merge \
	--depend=afterok:${JOB0##* }  \
	${output_dir}/slurm.bam_call_genotypes_merge)
echo $JOB1


