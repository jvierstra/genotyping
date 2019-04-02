#!/bin/bash

GZVCF_FILE=/net/seq/data/projects/genotyping/results.dgf-samples.merge2.genotype/filtered.all.hets-pass.recoded-final.dbSNP.vcf.gz
output_dir=/net/seq/data/projects/genotyping/results.dgf-samples.merge2.genotype/sample.counts

rm -rf ${output_dir}/logs && mkdir -p ${output_dir}/logs

module load bcftools/1.7

grep -v "\#" /net/seq/data/projects/genotyping/results.dgf-samples.merge2/sample.relatedness.cluster.txt \
| awk -v OFS="\t" '{ print $1, "/net/seq/data/projects/genotyping/results.dgf-samples.merge2/individual."$2".bam"; }' \
> ${output_dir}/inputs.txt

bcftools view -h ${GZVCF_FILE} | perl -ne 'print "$1\n" if /contig=<ID=([^,;]+),/' \
> ${output_dir}/chroms.txt

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
prefix=\`basename \${params[0]} | cut -d"." -f1\`

python2 /home/jvierstra/proj/code/genotyping/scripts/count_tags_bed.py ${GZVCF_FILE} \${params[0]} \${params[1]} f> ${output_dir}/\${prefix}.bed
__SCRIPT__

cat <<__SCRIPT__ > ${output_dir}/slurm.recode_vcf
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%A.%a.out
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=queue0

set -e -o pipefail

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

chrom=(\`cat ${output_dir}/chroms.txt | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1\`)

python2 /home/jvierstra/proj/code/genotyping/scripts/recode_vcf_with_tags_per_sample.py --chrom \${chrom} ${GZVCF_FILE} ${output_dir}/ ${output_dir}/inputs.txt  ${output_dir}/\${chrom}.vcf.gz

__SCRIPT__

cat <<__SCRIPT__ > ${output_dir}/slurm.merge
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%J.out
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=queue0
 
bcftools concat -Oz ${output_dir}/chr*.vcf.gz > ${output_dir}/merged.all.vcf.gz

#bcftools query -i 'GT="het"' -f '%CHROM\t%POS0\t%ID\t%REF\t%ALT\t[%ARD{0},]\t[%ARD{1},]\n' ${output_dir}/merged.all.vcf.gz \
#| python2 /home/jvierstra/proj/code/genotyping/scripts/compute_ai.py > ${output_dir}/merged.all.ai.bed

python /home/jvierstra/proj/code/genotyping/scripts/count_genotypes.py ${output_dir}/merged.all.vcf.gz > ${output_dir}/merged.all.ai.bed
__SCRIPT__



#njobs=$(wc -l < ${output_dir}/inputs.txt)
#echo "Number of datasets: $njobs"
#JOB0=$(sbatch --export=ALL --parsable \
#	--job-name=allelic.reads \
#	--array=1-${njobs} \
#	${output_dir}/slurm.bam_count_tags_chunk)
#echo $JOB0


njobs=$(wc -l < ${output_dir}/chroms.txt)
echo "Number of chromosomes: $njobs"
JOB1=$(sbatch --export=ALL --parsable \
	--job-name=allelic.reads.recode_vcf \
	--array=1-${njobs} \
	${output_dir}/slurm.recode_vcf)
echo $JOB1

JOB2=$(sbatch --export=ALL --parsable \
	--job-name=allelic.reads.merge \
	--depend=afterok:${JOB1} \
	${output_dir}/slurm.merge)
echo $JOB2

