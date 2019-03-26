# genotyping
Call genotypes directly from epigenomics datasets

# outline of typical pipeline
1) call_genotypes.sh            : detect significant SNPs over all BAM files
2) run_relatedness.sh           : determine relationships to group close relatives
3) cluster_relatedness.sh       : visualize how many close relations there are
4) merge_bams_by_individual.sh  : merge BAM files of close relatives
5) call_genotypes.sh            : detect significant SNPs over revised set of individuals after merging
6) counts_tags_by_individual.sh : creates a bed-matrix precursor for the next step
7) count_genotypes.caller.sh    : bed-matrix stats for het/homozygous information over all individuals

Step 3 can often be skipped if tissue culture history is known (edits/clones should be grouped together as one individual)

Step 7 has ref/alt allele information (col4), often taken from overlapping TOPMED info.  Alternative methods could be used.
  and scripts/create-refalt-alleles.tcsh has some comments about that.
  (In hindsight, this is not necessary (TOPMED) because we use the same reference genome.  To be updated later.  The current
   approach covers in excess of 93% of the significant variants called.)

The pipeline/example/ dir shows caller script examples.

