#!/bin/env python

import sys
import numpy as np
from scipy.stats import binom_test

def format_counts(a):
	return ':'.join(map(str, a))


for line in sys.stdin:

	fields = line.strip().split('\t')

	chrom = fields[0]
	start = int(fields[1])
	end = int(fields[2])

	(ref, alt) = fields[3].split("/")

	n_homozygous_ref = 0
	n_homozygous_alt = 0
	n_heterozygous = 0

	homozygous_ref = np.zeros(2, dtype = int)
	homozygous_alt = np.zeros(2, dtype = int)
	heterozygous = np.zeros(2, dtype = int)

	total = 0

	for ds in fields[4:]:

		(genotype, nr, na) = ds.split(":")
		counts = np.array([nr, na], dtype = np.int)

		if genotype == "./.":
			continue

		(a0, a1) = genotype.split("/")

		if a0 == a1:
			if a0 == ref:
				n_homozygous_ref += 1
				homozygous_ref += counts

			elif a0 == alt:
				n_homozygous_alt += 1
				homozygous_alt += counts
		else:
			n_heterozygous += 1
			heterozygous += counts

	total = n_homozygous_ref + n_homozygous_alt + n_heterozygous

	try:
		imbalance = float(heterozygous[0]) / float(heterozygous[0]+heterozygous[1]) 
	except ZeroDivisionError:
		imbalance = np.nan

	if heterozygous[0]+heterozygous[1] >= 50:
		p = binom_test(heterozygous[:2], p = 0.5)
	else:
		p = np.nan

	print "\t".join(map(str, [chrom, start, end, "%s/%s" % (ref, alt) , total, 
		n_homozygous_ref, format_counts(homozygous_ref), 
		n_homozygous_alt, format_counts(homozygous_alt), 
		n_heterozygous, format_counts(heterozygous), imbalance, p]))





