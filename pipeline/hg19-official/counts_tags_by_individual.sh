#!/bin/bash -x

usage () {
  echo -e "Usage: $0 <gzvcf_file> <output-dir>" >&2
  exit 2
}

if [[ $# != 2 ]]; then
   usage
fi


module load bcftools/1.7

GZVCF_FILE=$1
output_dir=$2

rm -rf ${output_dir}/logs && mkdir -p ${output_dir}/logs

#Get samples
bcftools query -l ${GZVCF_FILE} > $output_dir/inputs.txt

bcftools view -h ${GZVCF_FILE} | perl -ne 'print "$1\n" if /contig=<ID=([^,;]+),/' \
> ${output_dir}/chroms.txt

cat <<__SCRIPT__ > ${output_dir}/slurm.bam_count_tags_chunk
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%A.%a.out
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=queue0

module load python/2.7.11
module load numpy
module load htslib/1.7
module load pysam

set -x -e -o pipefail

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

params=(\`cat ${output_dir}/inputs.txt | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1\`)
prefix=\`basename \${params[0]} | cut -d"." -f1\`

python2 /net/seq/data/projects/altius-100-donors/hg19/scripts2.hg19/genotyping/scripts/count_tags_bed.py "${GZVCF_FILE}" "\${params[0]}" "\${params[0]}" > ${output_dir}/\${prefix}.bed
bgzip -c ${output_dir}/\$prefix.bed > ${output_dir}/\$prefix.bed.gz
tabix -p bed ${output_dir}/\$prefix.bed.gz
__SCRIPT__


cat <<__SCRIPT__ > ${output_dir}/slurm.recode_vcf
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%A.%a.out
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=queue0

set -x -e -o pipefail

module load python/2.7.11
module load gcc/5.3.0
module load numpy
module load atlas-lapack
module load scipy

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

chrom=(\`cat ${output_dir}/chroms.txt | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1\`)

python2 /net/seq/data/projects/altius-100-donors/scripts/genotyping/scripts/recode_vcf_with_tags.py --chrom "\${chrom}" "${GZVCF_FILE}" "${output_dir}/" "${output_dir}/\${chrom}.vcf.gz"

__SCRIPT__

cat <<__SCRIPT__ > ${output_dir}/slurm.merge
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%J.out
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=queue0

module load bcftools/1.7
module load python/2.7.11
module load numpy
module load gcc
module load atlas-lapack
module load scipy
module load htslib/1.7

bcftools concat -Oz ${output_dir}/chr*.vcf.gz > ${output_dir}/merged.all.vcf.gz

bcftools query -i 'GT="het"' -f '%CHROM\t%POS0\t%ID\t%REF\t%ALT\t[%ARD{0},]\t[%ARD{1},]\n' ${output_dir}/merged.all.vcf.gz \
| python2 /net/seq/data/projects/altius-100-donors/scripts/genotyping/scripts/compute_ai.py > ${output_dir}/merged.all.ai.bed

tabix -p vcf ${output_dir}/merged.all.vcf.gz
__SCRIPT__

njobs=$(wc -l < ${output_dir}/inputs.txt)
echo "Number of datasets: $njobs"
JOB0=$(sbatch --export=ALL --parsable \
	--job-name=allelic.reads \
	--array=1-${njobs} \
	${output_dir}/slurm.bam_count_tags_chunk)
echo $JOB0


njobs=$(wc -l < ${output_dir}/chroms.txt)
echo "Number of chromosomes: $njobs"
JOB1=$(sbatch --export=ALL --parsable \
	--job-name=allelic.reads.recode_vcf \
	--array=1-${njobs} \
	--depend=afterok:${JOB0} \
	${output_dir}/slurm.recode_vcf)
echo $JOB1


JOB2=$(sbatch --export=ALL --parsable \
	--job-name=allelic.reads.merge \
	--depend=afterok:${JOB1} \
	${output_dir}/slurm.merge)
echo $JOB2

