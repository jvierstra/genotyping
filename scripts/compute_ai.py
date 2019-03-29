#!/bin/env python

import sys
import numpy as np
from scipy.stats import binom_test

for line in sys.stdin:
	
	(chrom, pos, rsnum, ref, alt, refcounts,altcounts) = line.strip().split("\t")


	NR=np.array(map(int, refcounts.strip(',').split(',')), dtype=int)
	NA=np.array(map(int, altcounts.strip(',').split(',')), dtype=int)
	
	assert(len(NR) == len(NA))

	tref=np.sum(NR)
	talt=np.sum(NA)

	try:
		r = float(tref)/float(tref+talt)
	except ZeroDivisionError:
		r = np.nan

	if len(NR)<2 or (tref+talt)<50:
		p = np.nan
	else:
		p = binom_test([float(tref), float(talt)], p=0.5)


	print '\t'.join(map(str, [chrom, pos, int(pos)+1, rsnum, ref, alt, len(NR), tref, talt, r, p]))