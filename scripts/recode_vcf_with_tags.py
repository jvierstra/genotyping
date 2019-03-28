#!/bin/env python2

import sys
import logging
from os.path import basename, splitext, join
from argparse import ArgumentParser

import pysam
from scipy.stats import binom_test
import numpy as np

def parse_options(args):

	parser = ArgumentParser(description = "Update VCF file with filtered per-individual read counts")

	parser.add_argument("input_file", metavar = "input_file", type = str,
						help = "Path to VCF-format genotyping file. Format fields must include GT,GQ,PL,AD,DP tags")

	parser.add_argument("input_dir", metavar = "input_dir", type = str,
						help = "Input directory containing per-individual allelic tag counts")

	parser.add_argument("output_file", metavar = "output_file", type = str,
						help = "Output VCF file")

	parser.add_argument("--chrom", metavar = "chrom", type = str,
						default=None, help = "Restrict to a specific chromosome")


	return parser.parse_args(args)

class GenotypeError(Exception):
	pass

def main(argv = sys.argv[1:]):

	args = parse_options(argv)

	infile=pysam.VariantFile(args.input_file, mode='r', ignore_truncation=True)
	infile.header.formats.add("ARD", 'R', "Integer", "Allelic read depth filtered from BAM file")

	#out file with new VCF header INFO and FORMAT fields
	outfile=pysam.VariantFile(args.output_file, mode='w', header=infile.header)

	#outfile.header.info.add("GTC", 'G', "Integer", "Genotype counts")
	#outfile.header.info.add("AI", 4, "Float", "Total allelic reads (for heterozygotes), allelic ratio and binomial p-value")

	outfile.header.add_meta("ARD_command", "recode_vcf_with_tags.py " + ' '.join(["%s=%s" %(k, v) for k, v in vars(args).items()]))
	#outfile.close()
	#sys.exit()
	
	sample_infiles = dict()

	for sample in infile.header.samples:
		basefilepath = splitext(basename(sample))[0]
		sample_infiles[sample] = pysam.TabixFile(join(args.input_dir, basefilepath + ".bed.gz"))

	for var in infile.fetch(contig=args.chrom):

		outvar = var.copy()

		for sample, info in outvar.samples.items():
			
			try:
				if not all(info.alleles):
					raise StopIteration

				# try to read a line from the individual passed genotypes file
				# if we find a line, nothing changes
				# if no line exists that means we filtered the call
				tabix = sample_infiles[sample]
				row = tabix.fetch(var.contig, var.start, var.start+1, parser=pysam.asTuple()).next()
				
				#sanity check
				(a0, a1) = row[3].split('/')
				if a0!=info.alleles[0] or a1!=info.alleles[1]:
					raise GenotypeError()

				counts = row[5:7]
			
			except GenotypeError:
				logging.critical("Genotypes from VCF file do not match individual calls. Wrong file?")
				sys.exit(1)
			except StopIteration:
				counts = [0, 0]
	
			finally:
				info["ARD"] = map(int, counts)

		outfile.write(outvar)

	for sample, filehandle in sample_infiles.items():
		filehandle.close()

	outfile.close()
	infile.close()
		
if __name__ == "__main__":
	sys.exit(main())

