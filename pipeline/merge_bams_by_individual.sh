#!/bin/bash


cluster_file=/home/jvierstra/proj/genotypes/cluster.sample.bam.txt
output_dir=/home/jvierstra/proj/genotypes/results.individual

rm -rf ${output_dir}/logs && mkdir -p ${output_dir}/logs

njobs=$(cut -f2 ${cluster_file} | sort -g | uniq | tail -n 1)
echo ${njobs}

cat <<__SCRIPT__ > ${output_dir}/slurm.merge_bams
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%A.%a.out
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --partition=queue0

set -e -o pipefail

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

module load samtools/1.7

awk -v cl=\${SLURM_ARRAY_TASK_ID} '\$2==cl { print \$3; }' ${cluster_file} > \${TMPDIR}/bamfiles.txt; \

nfiles=\$(wc -l < \${TMPDIR}/bamfiles.txt)
first_file=\$(head -n 1 \${TMPDIR}/bamfiles.txt)

if [ \$nfiles -eq 1 ]; then
	cp \$first_file ${output_dir}/individual.\${SLURM_ARRAY_TASK_ID}.bam
else
	samtools merge -f -b \${TMPDIR}/bamfiles.txt -@2 ${output_dir}/individual.\${SLURM_ARRAY_TASK_ID}.bam
fi

samtools index ${output_dir}/individual.\${SLURM_ARRAY_TASK_ID}.bam

__SCRIPT__

JOB0=$(sbatch --export=ALL \
	--job-name=merging \
	--array=1-${njobs} \
	${output_dir}/slurm.merge_bams)
echo $JOB0