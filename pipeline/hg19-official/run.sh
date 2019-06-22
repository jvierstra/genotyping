#!/bin/bash -x

cd "$(dirname "$0")"

module purge

export script_dir=$(readlink -f genotyping/scripts/)

wait_for() {
  jobs=$("$@" | awk '/Submitted/ {print $NF}')
  while squeue | awk '!/Satisfied/ {print $1}' | grep -q -f <(echo $jobs) ; do
    echo "waiting for jobs: $jobs"
    sleep 60
  done
}

input_dir=$(readlink -f ../data/input.no36bpfilter)
output_dir=$(readlink -f ../data/output.no36bpfilter)
mkdir -p $output_dir/individual.counts

call_jobs=$(./call_genotypes.sh ${input_dir} ${output_dir})

# Wait for call_jobs to be complete
while squeue -u $(whoami) | awk '{print $1}' | grep -q -f <(echo "$call_jobs" | xargs -n1) ; do
  echo "waiting for jobs: $call_jobs"
  sleep 60
done

srun ./filter_hetzygs_sites_no36bp-filt.sh ${output_dir} #--dependency=afterok:$call_job

gzvcf=${output_dir}/filtered.all.hets-pass.recoded-final.vcf.gz
./counts_tags_by_individual.sh ${gzvcf} ${output_dir}/individual.counts

exit 0
