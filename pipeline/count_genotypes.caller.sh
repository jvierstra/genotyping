#!/bin/bash -x

usage () {
  echo -e "Usage: $0 <input-bed-mtx> <output>" >&2
  exit 2
}

if [[ $# != 2 ]]; then
   usage
fi

input=$1
output=$2

output_dir=$(dirname $output)

cat <<__SCRIPT__ > ${output_dir}/slurm.bam_count_genotypes
#!/bin/bash -x
#
#SBATCH --output=${output_dir}/logs/%A.%a.out
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=queue0

set -e -o pipefail

module add anaconda/2.1.0-dev
conda create -y -n conda-env"-$$" python=2.7
source activate conda-env"-$$"
conda install -y numpy
conda install -y scipy

python count_genotypes.py < ${input} > ${output}
__SCRIPT__

chmod 755 ${output_dir}/slurm.bam_count_genotypes

JOB0=$(sbatch --export=ALL \
    --job-name=genotype.counts \
    ${output_dir}/slurm.bam_count_genotypes)
echo $JOB0

exit 0
