#!/bin/env python2

import sys
import logging
from os.path import basename, splitext, join
from argparse import ArgumentParser

import pysam
from scipy.stats import binom_test
import numpy as np

def parse_options(args):

	parser = ArgumentParser(description = "Perform per sample filtering of a VCF file to get heterozygous sites")

	parser.add_argument("input_file", metavar = "input_file", type = str,
						help = "Path to VCF-format genotyping file. Format fields must include GT,GQ,PL,AD,DP tags")

	parser.add_argument("input_dir", metavar = "input_dir", type = str,
						help = "Output directory for per sample BED files with allelic tags counts")

	parser.add_argument("output_file", metavar = "output_file", type = str,
						help = "Output directory for per sample genotype BED files")

	parser.add_argument("--min_het_reads", metavar = "min_het_reads", type = int, default=50, 
						help = "Mininum number of HETEROZYGOUS reads to test for imbalance")

	parser.add_argument("--min_het_samples", metavar = "min_het_samples", type = int, default=1, 
						help = "Mininum number of HETEROZYGOUS samples to test for imbalance")

	parser.add_argument("--p", metavar = "p", type = float, default=0.5, 
						help = "Binomial paramter p used in imbalance test (default: %(default)s)")


	return parser.parse_args(args)

class GenotypeError(Exception):
	pass

def main(argv = sys.argv[1:]):

	args = parse_options(argv)

	print args

	infile=pysam.VariantFile(args.input_file, mode='r', ignore_truncation=True)
	infile.header.formats.add("ARD", 'R', "Integer", "Allelic read depth filtered from BAM file")

	#out file with new VCF header INFO and FORMAT fields
	outfile=pysam.VariantFile(args.output_file, mode='w', header=infile.header)

	#outfile.header.info.add("GTC", 'G', "Integer", "Genotype counts")
	#outfile.header.info.add("AI", 4, "Float", "Total allelic reads (for heterozygotes), allelic ratio and binomial p-value")

	outfile.header.add_meta("allelic_reads_command", "recode_vcf_with_tags.py " + ' '.join(["%s=%s" %(k, v) for k, v in vars(args).items()]))

	#outfile.close()
	#sys.exit()
	
	sample_infiles = dict()

	for sample in infile.header.samples:
		basefilepath = splitext(basename(sample))[0]
		sample_infiles[sample] = pysam.TabixFile(join(args.input_dir, basefilepath + ".bed.gz"))

	for var in infile.fetch():

		outvar = var.copy()

		hom_ref = np.zeros(2, dtype=int)
		hom_alt = np.zeros(2, dtype=int)
		het = np.zeros(2, dtype=int)

		n_hom_ref = n_hom_alt = n_het = 0

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

				counts = np.array(row[5:7], dtype=int)

				# set FORMAT TAG
			
				if a0==a1:
					if a0==var.ref:
						hom_ref += counts
						n_hom_ref += 1
					elif a0==var.alts[0]:
						hom_alt += counts
						n_hom_alt += 1
				else:
					het += counts
					n_het += 1
			
			except GenotypeError, e:
				logging.critical("Genotypes from VCF file do not match individual calls. Wrong file?")
				sys.exit(1)
			except StopIteration:
				counts = [0, 0]
	
			finally:
				info["ARD"] = map(int, counts)
		# compute allelic ratio
		#try:
		#	het_ratio = float(het[0]/float(np.sum(het)))
		#except ZeroDivisionError:
		#	het_ratio = -1

		#
		#if np.sum(het)>=args.min_het_reads and n_het>=args.min_het_samples:
		#	p = binom_test(het, p=args.p)
		#else:
		#	p = -1

		# SET INFO TAGS
		#outvar.info["AI"] = [het_ratio, -np.log(p)]
		#outvar.info["GTC"] = [n_hom_ref, n_hom_alt, n_het]

		outfile.write(outvar)

	for sample, filehandle in sample_infiles.items():
		filehandle.close()

	outfile.close()
	infile.close()
		
if __name__ == "__main__":
	sys.exit(main())

