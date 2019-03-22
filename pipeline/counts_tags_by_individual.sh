#!/bin/bash

module load bcftools/1.7

FASTA_CHROM_FILE=/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.chrom_sizes.bed
FASTA_FILE=/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.fa
GZVCF_FILE=/net/seq/data/projects/genotyping/results.dgf-samples.merge2.genotype/filtered.hwe.0.01.vcf.gz

output_dir=/net/seq/data/projects/genotyping/results.dgf-samples.merge2.genotype/individual.counts
rm -rf ${output_dir}/logs && mkdir -p ${output_dir}/logs

#Get samples
bcftools query -l ${GZVCF_FILE} > /tmp/samples.txt
cat /tmp/samples.txt > ${output_dir}/inputs.txt

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

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

params=(\`cat ${output_dir}/inputs.txt | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1\`)
prefix=\`basename \${params[0]} | cut -d"." -f1,2\`

python /home/jvierstra/proj/code/genotyping/scripts/count_tags.py ${GZVCF_FILE} \${params[0]} \${params[0]} > ${output_dir}/\${prefix}.bed
__SCRIPT__


JOB0=$(sbatch --export=ALL --parsable\
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

JOB1=$(sbatch --export=ALL --parsable\
	--job-name=allelic.reads.merge \
	--depend=afterok:${JOB0}  \
	${output_dir}/slurm.bam_count_tags_merge)
echo $JOB1
