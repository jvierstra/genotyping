#!/bin/env

FASTA_CHROM_FILE=/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.chrom_sizes.bed
FASTA_FILE=/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.fa
GZVCF_FILE=/home/jvierstra/proj/genotypes/results/filtered.all.bcftools.vcf.gz

output_dir=/home/jvierstra/proj/genotypes/results.allele-counts
rm -rf ${output_dir}/logs && mkdir -p ${output_dir}/logs


#Get samples
bcftools query -l ${GZVCF_FILE} > /tmp/samples.txt

cat /tmp/samples.txt > ${output_dir}/inputs.txt


# Make chunks
#cat ${FASTA_CHROM_FILE} \
#| grep -v random | grep -v _hap | grep -v chrM | grep -v chrUn | grep -v chrX | grep -v chrY | grep -v chrEBV \
#| cut -f1,3 \
#| awk -v step=150000000 -v OFS="\t" '{ for(i=step; i<=$2; i+=step) { print $1":"i-step+1"-"i; } print $1":"i-step+1"-"$2; }' \
#> /tmp/chunks.txt

#Make a samples x regions parameters file
#awk -v OFS="\t" 'FNR==NR{a[++n]=$0; next} {for(i=1; i<=n; i++) print $0, a[i]}' /tmp/samples.txt /tmp/chunks.txt > ${output_dir}/inputs.txt


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

#python /home/jvierstra/proj/genotypes/count_tags.py ${GZVCF_FILE} \${params[1]} \${params[0]} > ${output_dir}/\${prefix}.\${params[0]}.bed
python /home/jvierstra/proj/genotypes/count_tags.py ${GZVCF_FILE} \${params[0]} > ${output_dir}/\${prefix}.bed

__SCRIPT__


JOB0=$(sbatch --export=ALL \
	--job-name=allelic.reads \
	--array=1-${njobs} \
	${output_dir}/slurm.bam_count_tags_chunk)
echo $JOB0

cmd="";\
LNS=(`cat inputs.txt | xargs -L1 'basename' | cut -f1 -d"."`);\
for ln in ${LNS[@]}; do cmd="$cmd <(cut -f5 ${ln}.bed)"; done;\
cmd="paste <(cut -f1-4 ${LNS[0]}.bed) $cmd"; eval $cmd > merged.counts.no-filter.bed

