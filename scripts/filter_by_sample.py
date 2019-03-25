#!/bin/env python2

import sys
import logging
from os.path import basename, splitext, join
from argparse import ArgumentParser

import pysam


def parse_options(args):

	parser = ArgumentParser(description = "Perform per sample filtering of a VCF file to get heterozygous sites")

	parser.add_argument("vcf_file", metavar = "vcf_file", type = str,
						help = "Path to VCF-format genotyping file. Format fields must include GT,GQ,PL,AD,DP tags.")

	parser.add_argument("output_dir", metavar = "output_dir", type = str,
						help = "Output directory for per sample genotype BED files")

	parser.add_argument("--min_dp", metavar = "min_dp", type = int, default=12, 
						help = "Mininum read depth over site (default: %(default)s)")

	parser.add_argument("--min_ad", metavar = "min_ad", type = int, default=4, 
						help = "Mininum reads per allele over heterozygous sites (default: %(default)s)")

	parser.add_argument("--min_gq", metavar = "min_gq", type = int, default=50, 
						help = "Mininum genotype quality score (default: %(default)s)")

	return parser.parse_args(args)


def main(argv = sys.argv[1:]):

	args = parse_options(argv)

	vcf=pysam.VariantFile(args.vcf_file)

	sample_outfiles = dict()

	for sample in vcf.header.samples:
		basefilepath = splitext(basename(sample))[0]
		sample_outfiles[sample] = open(join(args.output_dir, basefilepath + ".bed"), 'w')

	gtmap=[None, None]

	for var in vcf:

		n_hets=0
		n_homs=0

		gtmap[0]=var.ref
		gtmap[1]=var.alts[0]

		for sample, info in var.samples.items():

			# if no gt called omit
			gt=info["GT"]
			if not any(gt):
				continue
			
			ref=gtmap[gt[0]]
			alt=gtmap[gt[1]]

			#filter on total genotype depth
			dp=info["DP"]
			if dp<args.min_dp:
				continue

			# filter on allele depth also (only for hets)
			ad=info["AD"]
			if ref!=alt and (ad[0]<args.min_ad or ad[1]<args.min_ad):
				continue

			# genotype quality
			gq=info["GQ"]
			if gq<args.min_gq:
				continue

			outline="\t".join(map(str, [var.contig, var.start, var.start+1, ref, alt, gq, dp, str(ad[0]) + ":" + str(ad[1])]))
			sample_outfiles[sample].write(outline + "\n")

			if ref==alt:
				n_homs+=1
			else:
				n_hets+=1

		print "\t".join(map(str, [var.contig, var.start, var.start+1, var.ref, var.alts[0], var.qual, var.info["DP"], n_homs, n_hets]))


	for sample, filehandle in sample_outfiles.items():
		filehandle.close()


    
if __name__ == "__main__":
    sys.exit(main())

