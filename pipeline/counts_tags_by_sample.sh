#!/bin/bash

GZVCF_FILE=/net/seq/data/projects/genotyping/results.dgf-samples.merge2.genotype/filtered.all.hets-pass.recoded-final.vcf.gz
output_dir=/net/seq/data/projects/genotyping/results.dgf-samples.merge2.genotype/sample.counts

grep -v "\#" /net/seq/data/projects/genotyping/results.dgf-samples.merge2/sample.relatedness.cluster.txt \
| awk -v OFS="\t" '{ print $1, "/net/seq/data/projects/genotyping/results.dgf-samples.merge2/individual."$2".bam"; }' \
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

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

params=(\`cat ${output_dir}/inputs.txt | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1\`)
prefix=\`basename \${params[0]} | cut -d"." -f1\`

python /home/jvierstra/proj/code/genotyping/scripts/count_tags_bed.py ${GZVCF_FILE} \${params[0]} \${params[1]} > ${output_dir}/\${prefix}.bed
__SCRIPT__

cat <<__SCRIPT__ > ${output_dir}/slurm.bam_count_tags_merge
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%J.out
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=queue0

module load bcftools/1.7

cmd=""
LNS=(\`cut -f1 ${output_dir}/inputs.txt | xargs -L1 'basename' | cut -f1 -d"."\`)
for ln in \${LNS[@]}; do 
	cmd="\$cmd <(cut -f4,6,7 \${ln}.bed | tr '\t' ':')"
done;
#cmd="paste <(cut -f1-3 \${LNS[0]}.bed) $cmd"
cmd="paste <(bcftools query -f '%CHROM\t%POS0\t%POS\t%REF/%ALT\n' ${GZVCF_FILE}) \$cmd"
eval \$cmd > ${output_dir}/merged.counts.no-filter.bed
__SCRIPT__

rm -rf ${output_dir}/logs && mkdir -p ${output_dir}/logs
JOB0=$(sbatch --export=ALL --parsable \
	--job-name=allelic.reads \
	--array=1-${njobs} \
	${output_dir}/slurm.bam_count_tags_chunk)
echo $JOB0

#JOB1=$(sbatch --export=ALL --parsable\
#	--job-name=allelic.reads.merge \
#	--depend=afterok:${JOB0}  \
#	${output_dir}/slurm.bam_count_tags_merge)
#echo $JOB1
