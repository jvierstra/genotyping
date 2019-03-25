#!/bin/env python2

import sys
import logging
from os.path import basename, splitext, join
from argparse import ArgumentParser

import pysam


def parse_options(args):

	parser = ArgumentParser(description = "Perform per sample filtering of a VCF file to get heterozygous sites")

	parser.add_argument("input_file", metavar = "input_file", type = str,
						help = "Path to VCF-format genotyping file. Format fields must include GT,GQ,PL,AD,DP tags.")

	parser.add_argument("input_dir", metavar = "input_dir", type = str,
						help = "Output directory for per sample genotype BED files")

	parser.add_argument("output_file", metavar = "output_file", type = str,
						help = "Output directory for per sample genotype BED files")

	return parser.parse_args(args)


def main(argv = sys.argv[1:]):

	args = parse_options(argv)

	infile=pysam.VariantFile(args.input_file, ignore_truncation=True)

	outfile=pysam.VariantFile(args.output_file, mode='w', header=infile.header)

	sample_infiles = dict()

	for sample in infile.header.samples:
		basefilepath = splitext(basename(sample))[0]
		sample_infiles[sample] = pysam.TabixFile(join(args.input_dir, basefilepath + ".bed.gz"))

	for var in infile.fetch():

		outvar = var.copy()

		for sample, info in outvar.samples.items():
			
			# never had a genotype called
			if not all(info.alleles):
				continue

			# try to read a line from the individual passed genotypes file
			# if we find a lin, nothing changes
			# if no line exists that means we filtered the call
			tabix = sample_infiles[sample]
			try:
				row = tabix.fetch(var.contig, var.start, var.start+1, parser=pysam.asTuple()).next()
			except StopIteration:
				info.alleles = [None, None]

		outfile.write(outvar)

	for sample, filehandle in sample_infiles.items():
		filehandle.close()

	outfile.close()
	infile.close()
	    
if __name__ == "__main__":
    sys.exit(main())

