#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

# download and compile phaser from https://github.com/secastel/phaser.git
# download and gunzip hg38 files required by phaser
#  https://www.dropbox.com/s/1zapu5n4aeyi1g6/hg38_hla.chr.bed.gz?dl=0
#  https://www.dropbox.com/s/9cn9477bcutvuc7/hg38_haplo_count_blacklist.chr.bed.gz?dl=0

set dstamp = 2019-04-04
set script = $HOME/Github/phaser/phaser/phaser.py
set bamd = ../data/$dstamp/bams
set donorsmatched = ../data/$dstamp/prodQC/donors.matched.simple # dnase-rnase matches
set haplos = ../data/hg38_hla.chr.bed
set haplo_blacklist = ../data/hg38_haplo_count_blacklist.chr.bed
set vcfd = ../results/$dstamp/dnaseI/output/phasing
set vcfds = (`find $vcfd -maxdepth 1 -mindepth 1 -type d`)
set baseoutd = ../results/$dstamp/rnase/output/imputed

source /net/module/Modules/default/tcsh

module add bcftools/1.7
module add anaconda/2.1.0-dev

set envnm = conda-env-impute-rnase-t3
set activator = `which conda`
set activator = $activator:h/activate

@ exists = `conda env list | awk -v en=$envnm '$1 == en' | wc -l`
if ( $exists == 0 ) then
/bin/bash -x <<__EOF__
  module add anaconda/2.1.0-dev
  conda create -y -n $envnm python=2.7

  conda install -y -n $envnm libgfortran
  conda install -y -n $envnm scipy
  conda install -y -n $envnm pysam

  # for phaser_gene.ae.py later on
  conda install -y -n $envnm pandas
  conda install -c conda-forge -y -n $envnm intervaltree 
__EOF__
endif

foreach d ($vcfds)
  set vcfnm = $d:t

  foreach chr (`find $d/ -maxdepth 1 -mindepth 1 -type f -name "*.vcf.gz"`)
    foreach sample (`bcftools query -l $chr`)
      set samplenm = `echo $sample | tr '/' '.' | cut -f2- -d'.'`
      set dnasenm = $sample:t:r
      set rnasenm = `awk -v d=$dnasenm '$3 == d { print $4 }' $donorsmatched`
      if ( $rnasenm == "" ) then
        printf "No-match %s\n" $dnasenm
        continue
      endif

      set bam = (`find -L $bamd/ -mindepth 1 -maxdepth 1 -type f -name "*-$rnasenm.bam"`)
      if ( $#bam != 1 ) then
        printf "Problem with %s\n" $rnasenm
        exit -1
      endif

      set chrnm = `echo $chr:t | cut -f1 -d'.'`
      set output_dir = $baseoutd/$vcfnm/$samplenm
      mkdir -p $output_dir/logs

      sbatch --job-name=phaser.$vcfnm.$samplenm <<__SCRIPT__
#!/bin/bash -x
#
#SBATCH --output=$output_dir/logs/$chrnm.%A.%a.out
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=queue1

        module add anaconda/2.1.0-dev
        source $activator $envnm

        module add samtools/1.7
        module add htslib/1.7
        module add bedtools
        module add bcftools/1.7

        python $script \
                --vcf $chr \
                --bam $bam \
                --paired_end 1 \
                --mapq 255 \
                --baseq 10 \
                --sample $sample \
                --blacklist $haplos \
                --haplo_count_blacklist $haplo_blacklist \
                --threads 4 \
                --o $output_dir/$chrnm
__SCRIPT__
    end
  end
end

exit 0
