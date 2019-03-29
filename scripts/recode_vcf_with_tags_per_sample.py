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

	parser.add_argument("variant_file", metavar = "variant_file", type = str,
						help = "Path to VCF-format genotyping file.")

	parser.add_argument("input_dir", metavar = "input_dir", type = str,
						help = "Input directory containing per-sample allelic tag counts")

	parser.add_argument("sample_map_file", metavar = "sample_map_file", type = str,
						help = "Sample to individual mapping file")

	parser.add_argument("output_file", metavar = "output_file", type = str,
						help = "Output VCF file")

	parser.add_argument("--chrom", metavar = "chrom", type = str,
					default=None, help = "Restrict to a specific chromosome")


	return parser.parse_args(args)

class GenotypeError(Exception):
	pass

def main(argv = sys.argv[1:]):

	args = parse_options(argv)

	samples={}
	with open(args.sample_map_file) as f:
		for line in f:
			(fname, individual) = line.strip().split("\t")
			dataset=basename(fname).split('.')[0]
			samples[dataset] = individual

	#print samples

	infile=pysam.VariantFile(args.variant_file, mode='r', ignore_truncation=True)
	

	#out file with new VCF header INFO and FORMAT fields
	header=pysam.VariantHeader()
	for record in infile.header.records:
		header.add_record(record)
	
	header.formats.clear_header()
	header=header.copy() # hack because the clearing the formats doesn't actually delete them "behind the scenes"

	header.formats.add('GT', 1, "Integer", "Genotype")
	header.formats.add('ARD', 'R', 'Integer', 'Allelic read depth')

	#
	sample_infiles = dict()
	for samp in samples.keys():
		header.add_sample(samp)

		sample_infiles[samp] = pysam.TabixFile(join(args.input_dir, samp + ".bed.gz"))

	outfile=pysam.VariantFile(args.output_file, mode='w', header=header)
	outfile.header.add_meta("ARD_command", "recode_vcf_with_tags_per_sample.py " + ' '.join(["%s=%s" %(k, v) for k, v in vars(args).items()]))

	#
	for var in infile.fetch(contig=args.chrom):

		hom_ref = np.zeros(2, dtype=int)
		hom_alt = np.zeros(2, dtype=int)
		het = np.zeros(2, dtype=int)

		n_hom_ref = n_hom_alt = n_het = 0
		
		samples_format = []

		for samp, indiv in samples.items():

			alleles = var.samples[indiv].alleles
			gt = var.samples[indiv]['GT']

			try:
				
				if not all(alleles):
					raise StopIteration

				# try to read a line from the individual passed genotypes file
				# if we find a line, nothing changes
				# if no line exists that means we filtered the call
				tabix = sample_infiles[samp]
				row = tabix.fetch(var.contig, var.start, var.start+1, parser=pysam.asTuple()).next()
				
				#sanity check
				(a0, a1) = row[3].split('/')
				if a0!=alleles[0] or a1!=alleles[1]:
					raise GenotypeError()

				counts = row[5:7]

			except GenotypeError:
				logging.critical("Genotypes from VCF file do not match individual calls. Wrong file?")
				sys.exit(1)
			except StopIteration:
				counts = [0, 0]
			finally:
				samples_format.append({ 'GT': gt,  'ARD': map(int,counts)})

		outvar = outfile.new_record(contig=var.contig, start=var.start, stop=var.stop, alleles=var.alleles, id=var.id, qual=var.qual, filter=var.filter, info=var.info, samples=samples_format)
		outfile.write(outvar)

	for sample, filehandle in sample_infiles.items():
		filehandle.close()

	outfile.close()
	infile.close()
		
if __name__ == "__main__":
	sys.exit(main())

