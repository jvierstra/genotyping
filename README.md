# genotyping
Call genotypes directly from epigenomics datasets

# outline of typical pipeline
1) call_genotypes.sh           : detect significant SNPs over all BAM files
2) run_relatedness.sh          : determine relationships to group close relatives
3) cluster_relatedness.sh      : visualize how many close relations there are
4) merge_bams_by_individual.sh : merge BAM files of close relatives
5) call_genotypes.sh           : detect significant SNPs over revised set of individuals after merging

The pipeline/example/ dir shows caller script examples.

