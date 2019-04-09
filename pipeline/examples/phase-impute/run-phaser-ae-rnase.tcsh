#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-04
set script = ~sjn/Github/phaser/phaser_gene_ae/phaser_gene_ae.py
set genes = ../data/genes.bed
set baseind = ../results/$dstamp/rnaseq/output/imputed
set baseoutd = ../results/$dstamp/rnaseq-ae/output/expression

source /net/module/Modules/default/tcsh

module add anaconda/2.1.0-dev

set envnm = conda-env-impute-rnaseseq-t3
set activator = `which conda`
set activator = $activator:h/activate

foreach vcf (`find $baseind/ -maxdepth 1 -mindepth 1 -type d`)
  set vcfnm = $vcf:t
  foreach sample (`find $vcf/ -maxdepth 1 -mindepth 1 -type d`)
    set samplenm = $sample:t
    foreach counts (`find $sample/ -maxdepth 1 -mindepth 1 -type f -name "*.haplotypic_counts.txt"`)
      set chrnm = `echo $counts:t | cut -f1 -d'.'`
      set output_dir = $baseoutd/$vcfnm/$samplenm
      mkdir -p $output_dir/logs

      sbatch --job-name=phaser.$vcfnm.$samplenm <<__SCRIPT__
#!/bin/bash -x
#
#SBATCH --output=$output_dir/logs/$chrnm.%A.%a.out
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=queue1

      module add anaconda/2.1.0-dev
      source $activator $envnm

      module add samtools/1.7
      module add htslib/1.7
      module add bedtools
      module add bcftools/1.7

      python $script \
              --haplotypic_counts $counts \
              --features $genes \
              --o $output_dir/$chrnm.ae.fullannotation.txt

      awk -F"\t" '\$10 != ""' $output_dir/$chrnm.ae.fullannotation.txt \
       > $output_dir/$chrnm.ae.txt

      rm -f $output_dir/$chrnm.ae.fullannotation.txt
__SCRIPT__
    end
  end
end

exit 0
