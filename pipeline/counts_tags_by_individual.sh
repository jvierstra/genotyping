#!/bin/bash

usage () {
  echo -e "Usage: $0 <genome-chromInfo.bed> <genome.fa> <hwe.0.01.vcf.gz> <output-dir>" >&2
  exit 2
}

if [[ $# != 4 ]]; then
   usage
fi

FASTA_CHROM_FILE=$1
FASTA_FILE=$2
GZVCF_FILE=$3
output_dir=$4

module load python
module load bcftools

#Get samples
bcftools query -l ${GZVCF_FILE} \
  | awk '{ n=split($1, a, "/"); nm=a[7]"."a[n]; print $1"\t"nm; }' \
  > ${output_dir}/inputs.txt

njobs=$(wc -l < ${output_dir}/inputs.txt)
echo "njobs: $njobs"

cat <<__SCRIPT__ > ${output_dir}/slurm.bam_count_tags_chunk
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%A.%a.out
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=queue0

set -e -o pipefail

module add numpy

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

params=(\`cat ${output_dir}/inputs.txt | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1 | cut -f1\`)
echo \$params

prefix=(\`cat ${output_dir}/inputs.txt | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1 | cut -f2\`)
echo \$prefix

python count_tags.py ${GZVCF_FILE} \${params[0]} \${params[0]} > ${output_dir}/\${prefix}.bed
__SCRIPT__


JOB0=$(sbatch --export=ALL \
	--job-name=allelic.reads \
	--array=1-${njobs} \
	${output_dir}/slurm.bam_count_tags_chunk)
echo $JOB0

cat <<__SCRIPT__ > ${output_dir}/slurm.bam_count_tags_merge
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%J.out
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=queue0

cmd=""
LNS=(\`cat ${output_dir}/inputs.txt | xargs -L1 'basename' | cut -f1,2 -d"."\`)
for ln in \${LNS[@]}; do 
	cmd="\$cmd <(cut -f5 \${ln}.bed)"
done;
cmd="paste <(cut -f1-4 \${LNS[0]}.bed) $cmd"
eval \$cmd > ${output_dir}/merged.counts.no-filter.bed
__SCRIPT__

JOB1=$(sbatch --export=ALL \
	--job-name=allelic.reads.merge \
	--depend=afterok:${JOB0##* }  \
	${output_dir}/slurm.bam_count_tags_merge)
echo $JOB1
