#!/bin/bash -x

usage () {
  echo -e "Usage: $0 <sigvars-vcf.gz> <output-dir>" >&2
  exit 2
}

if [[ $# != 2 ]]; then
   usage
fi

source /net/module/Modules/default/bash
module add bcftools/1.7

vcfgzfile=$1
output_dir=$2

mkdir -p ${output_dir}/logs

bcftools view -h ${vcfgzfile} | grep "contig" | awk '$1 !~ /_/' | perl -ne 'print "$1\n" if /contig=<ID=(chr[0-9XY]+)/' > ${output_dir}/chroms.txt
awk '{ print substr($1, 4), $1; }' ${output_dir}/chroms.txt > ${output_dir}/chroms.rename.txt


#bcftools view -h ALL.chr1_GRCh38.genotypes.20170504.bcf | grep "contig" | perl -ne 'print "$1\n" if /length=([0-9]+)/' > /tmp/chrom.lens.txt
#paste /tmp/chrom.names.txt /tmp/chrom.lens.txt | wc -l
#paste /tmp/chrom.names.txt /tmp/chrom.lens.txt \
#	| awk -v OFS="\t" -v winsize=20000000 \
#	'{ for(i=0; i<$2; i+=winsize) { end=(i+winsize)>$2?$2:i+winsize;  print $1, i, end; } }' \
#> {$output_dir}/chunks.txt


cat <<__SCRIPT__ > ${output_dir}/slurm.chunk
#!/bin/bash -x
#
#SBATCH --output=${output_dir}/logs/%A.%a.out
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=queue1

set -e -o pipefail

source /net/module/Modules/default/bash
module add bcftools/1.7
module add htslib/1.7
module add vcftools

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

chrom=(\`cat ${output_dir}/chroms.txt | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1\`)


## prephasing
/home/jvierstra/.local/src/Eagle_v2.4.1/eagle \
	--numThreads 16 \
	--outPrefix \${TMPDIR}/prephased \
	--geneticMapFile /home/jvierstra/.local/src/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
	--vcfRef /home/jvierstra/data/1k_genomes/refpanel/bcf/ALL.\${chrom}_GRCh38.genotypes.20170504.bcf \
	--vcfTarget ${vcfgzfile} \
	--vcfOutFormat u \
	--allowRefAltSwap \
	--noImpMissing \
	--chrom \${chrom}

(bcftools view --no-version -h \${TMPDIR}/prephased.bcf \
 | sed 's/^##contig=<ID=chr/##contig=<ID=/'; \
	bcftools view --no-version -H -c 2 \${TMPDIR}/prephased.bcf | sed 's/^chr//') \
 | bcftools view -Ov > \${TMPDIR}/prephased.chrom-fix.vcf

# produces unsorted vcf file (tabix complains) without 'chr' in contig names
/home/jvierstra/.local/src/Minimac3/bin/Minimac3-omp \
	--refHaps /home/jvierstra/data/1k_genomes/refpanel/m3vcfs/ALL.\${chrom}_GRCh38.genotypes.20170504.m3vcf.gz \
	--haps \${TMPDIR}/prephased.chrom-fix.vcf \
	--chr \${chrom#chr} \
	--prefix ${output_dir}/\${chrom}.imputed \
	--allTypedSites

# sort and add 'chr' to contig names; put samples/columns in a deterministic order
vcf-sort -p 4 ${output_dir}/\${chrom}.imputed.dose.vcf.gz \
 | bcftools annotate --rename-chrs ${output_dir}/chroms.rename.txt \
 | vcf-shuffle-cols -t ${vcfgzfile} \
 | bgzip -c > ${output_dir}/\${chrom}.imputed.dose.vcf.gz.new

mv ${output_dir}/\${chrom}.imputed.dose.vcf.gz.new ${output_dir}/\${chrom}.imputed.dose.vcf.gz
tabix -p vcf ${output_dir}/\${chrom}.imputed.dose.vcf.gz

__SCRIPT__

sbatch --export=ALL --parsable \
	--job-name=impute.genotype \
	--array=1-$(wc -l < ${output_dir}/chroms.txt) \
	 ${output_dir}/slurm.chunk

exit 0
