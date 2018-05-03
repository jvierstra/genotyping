#!/bin/env

###
# TODO:
#	-Check relatedness (vcf-tools)
#	-Filter SNPs
#

FASTA_CHROM_FILE=/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.chrom_sizes.bed
FASTA_FILE=/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.fa

output_dir=/net/seq/data/projects/genotyping/results.allsamples

rm -rf ${output_dir}/logs && mkdir -p ${output_dir}/logs

# Make chunks
cat ${FASTA_CHROM_FILE} \
| grep -v random | grep -v _hap | grep -v chrM | grep -v chrUn | grep -v chrX | grep -v chrY | grep -v chrEBV \
| cut -f1,3 \
| awk -v step=5000000 -v OFS="\t" '{ for(i=step; i<=$2; i+=step) { print $1":"i-step+1"-"i; } print $1":"i-step+1"-"$2; }' \
> ${output_dir}/inputs.txt

njobs=$(wc -l < ${output_dir}/inputs.txt)



# Make BAM file list
cut -f3  /home/jvierstra/proj/ftd.updated/hg38.DGF.ver4.params.txt > ${output_dir}/filelist.txt


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

bcftools mpileup -r \${region} -Q 20 -d 1000 -I -E -f ${FASTA_FILE} -b ${output_dir}/filelist.txt -Ou \
| bcftools call -f GQ -cv -Ov \
| vcftools --stdout \
	--minQ 500 --minGQ 30 --minDP 30 --max-alleles 2 \
	--recode --recode-INFO-all --vcf - \
| bgzip -c > ${output_dir}/\${region}.filtered.vcf.gz

#tabix -p vcf ${output_dir}/\${region}.filtered.vcf.gz

rm -rf \${TMPDIR}
__SCRIPT__

cat <<__SCRIPT__ > ${output_dir}/slurm.bam_call_genotypes_merge
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%A.%a.out
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=queue0

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

module load bcftools/1.7
module load vcftools/0.1.14
module load htslib/1.7

ls ${output_dir}/*.filtered.vcf.gz > \${TMPDIR}/files.txt

# Merge chunks
bcftools concat -f \${TMPDIR}/files.txt -Oz -o ${output_dir}/filtered.all.vcf.gz
tabix -p vcf ${output_dir}/filtered.all.vcf.gz

#
vcftools --temp \${TMPDIR} --gzvcf filtered.all.vcf.gz --out ${output_dir}/filtered.all --relatedness --het --depth --TsTv-by-count --TsTv-by-qual 

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


