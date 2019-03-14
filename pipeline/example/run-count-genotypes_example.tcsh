#!/bin/tcsh -efx

set pyscript = ../../scripts/count_genotypes.py
set caller = ../call_genotypes.sh
set indir = ../results.analysis.merged/results.genotype.inmtx
set inputs = (`find $indir/ -maxdepth 1 -mindepth 1 -type f -name "*.bed.mtx"`)

set output_dir = ../results.analysis.merged/results.genotype-counts
mkdir -p $output_dir/logs

cp $pyscript .
cp $caller .

foreach input ($inputs)
  set output = $output_dir/$input:t:r
  $caller $input $output
end

exit 0
