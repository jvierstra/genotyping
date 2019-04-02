#!/bin/bash

vcfgzfile=/net/seq/data/projects/altius-100-donors/archives/2019-03-29/output/filtered.all.vcf.gz
output_dir="/net/seq/data/projects/altius-100-donors/archives/2019-03-29/output/imputed"

mkdir -p ${output_dir}/logs

bcftools view -h ${vcfgzfile} | grep "contig" | perl -ne 'print "$1\n" if /contig=<ID=(chr[0-9XY]+)/' > ${output_dir}/chroms.txt


#bcftools view -h ALL.chr1_GRCh38.genotypes.20170504.bcf | grep "contig" | perl -ne 'print "$1\n" if /length=([0-9]+)/' > /tmp/chrom.lens.txt
#paste /tmp/chrom.names.txt /tmp/chrom.lens.txt | wc -l
#paste /tmp/chrom.names.txt /tmp/chrom.lens.txt \
#	| awk -v OFS="\t" -v winsize=20000000 \
#	'{ for(i=0; i<$2; i+=winsize) { end=(i+winsize)>$2?$2:i+winsize;  print $1, i, end; } }' \
#> {$output_dir}/chunks.txt



cat <<__SCRIPT__ > ${output_dir}/slurm.chunk
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%A.%a.out
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=queue0

set -e -o pipefail

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

chrom=(\`cat ${output_dir}/chroms.txt | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1\`)

## prepahseing
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


/home/jvierstra/.local/src/Minimac3/bin/Minimac3-omp \
	--refHaps /home/jvierstra/data/1k_genomes/refpanel/m3vcfs/ALL.\${chrom}_GRCh38.genotypes.20170504.m3vcf.gz \
	--haps \${TMPDIR}/prephased.chrom-fix.vcf \
	--chr \${chrom#chr} \
	--prefix ${output_dir}/\${chrom}.imputed \
	--allTypedSites

__SCRIPT__

sbatch --export=ALL --parsable \
	--job-name=impute.genotype \
	--array=1-$(wc -l < ${output_dir}/chroms.txt) \
	 ${output_dir}/slurm.chunk