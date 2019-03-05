#!/bin/bash -x

###
# TODO:
#	-Check relatedness (vcf-tools)
#	-Filter SNPs
#


usage () {
  echo -e "Usage: $0 <dir-of-bams> <output-dir> <chromInfo.bed> <genome.fasta>" >&2
  exit 2
}

if [[ $# != 4 ]]; then
   usage
fi

data_dir=$1
output_dir=$2
FASTA_CHROM_FILE=$3
FASTA_FILE=$4

rm -rf ${output_dir}/logs && mkdir -p ${output_dir}/logs

# Make chunks
cat ${FASTA_CHROM_FILE} \
| grep -v random | grep -v _hap | grep -v chrM | grep -v chrUn | grep -v chrX | grep -v chrY | grep -v chrEBV \
| cut -f1,3 \
| awk -v step=5000000 -v OFS="\t" '{ for(i=step; i<=$2; i+=step) { print $1":"i-step+1"-"i; } print $1":"i-step+1"-"$2; }' \
> ${output_dir}/inputs.txt

njobs=$(wc -l < ${output_dir}/inputs.txt)



# Make BAM file list
readlink `find -L $data_dir/ -maxdepth 1 -mindepth 1 -type f -name "*.bam"` > ${output_dir}/filelist.txt

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
	--minQ 500 --minGQ 50 --minDP 12 --max-alleles 2 \
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

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

module load bcftools/1.7
module load vcftools/0.1.14
module load htslib/1.7

cat ${output_dir}/inputs.txt | awk -v dir="${output_dir}" '{ print dir"/"\$1".filtered.vcf.gz"; }' > \${TMPDIR}/mergelist.txt
bcftools concat -f \${TMPDIR}/mergelist.txt -Oz -o ${output_dir}/filtered.all.vcf.gz
tabix -p vcf ${output_dir}/filtered.all.vcf.gz

vcftools --gzvcf ${output_dir}/filtered.all.vcf.gz --hwe 0.01 --stdout --recode --recode-INFO-all | bgzip -@2 -c > ${output_dir}/filtered.hwe.0.01.vcf.gz
tabix -p vcf ${output_dir}/filtered.hwe.0.01.vcf.gz

#bcftools query -f '%CHROM\t%POS0\t%POS\t%REF/%ALT\n' filtered.hwe.0.01.vcf.gz > filtered.hwe.0.01.bed
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

exit 0
