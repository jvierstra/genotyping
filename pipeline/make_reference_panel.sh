#!/bin/bash

#Make an imputation reference panel from 1000 genomes project phhase 3 release

output_dir="/home/jvierstra/data/1k_genomes/refpanel"
mkdir -p ${output_dir}/bcf ${output_dir}/m3vcfs ${output_dir}/logs

ls ~/data/1k_genomes/*.vcf.gz | perl -ne 'print "$1\n" if /(chr[0-9XY]+)/' > ${output_dir}/chroms.txt


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

module load bcftools/1.7

base="ALL.\${chrom}_GRCh38.genotypes.20170504"
vcfgzfile=/home/jvierstra/data/1k_genomes/\${base}.vcf.gz

tmp_vcffile=\${TMPDIR}/panel.vcf
tmp_vcffile_chromfix=\${TMPDIR}/panel.chromfix.vcf

panel_bcffile=${output_dir}/bcf/\${base}.bcf
panel_m3vcfsfile=${output_dir}/m3vcfs/\${base}

(bcftools view --no-version -h \${vcfgzfile} \
	| grep -v "^##contig=<ID=[GNh]" | sed 's/^##contig=<ID=MT/##contig=<ID=chrM/;s/^##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/'; \
	bcftools view --no-version -H -c 2 \${vcfgzfile} | sed 's/^/chr/') \
| bcftools norm --no-version -Ou -m -any \
| bcftools norm --no-version -Ou -d none -f /home/jvierstra/data/1k_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
| bcftools view -g^miss -Ov -o \${tmp_vcffile}


# Make BCF file for EAGLE2
bcftools view \${tmp_vcffile} -Ob -o \${panel_bcffile} && \
bcftools index -f \${panel_bcffile}

#Remove "chr" from contig line and compress to vcf.gz for minimac3
(bcftools view --no-version -h \${tmp_vcffile} \
	| sed 's/^##contig=<ID=chr/##contig=<ID=/'; \
	bcftools view --no-version -H -c 2 \${tmp_vcffile} | sed 's/^chr//') \
| bcftools view -Ov > \${tmp_vcffile_chromfix}

/home/jvierstra/.local/src/Minimac3/bin/Minimac3-omp \
	--refHaps \${tmp_vcffile_chromfix} \
	--processReference \
	--prefix \${panel_m3vcfsfile}

__SCRIPT__

sbatch --export=ALL --parsable \
	--job-name=make.ref.panel \
	--array=1-$(wc -l < ${output_dir}/chroms.txt) \
	 ${output_dir}/slurm.chunk
