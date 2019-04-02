import sys
import pysam
import numpy as np
from scipy.stats import binom_test

vcf_file=pysam.VariantFile(sys.argv[1])

for var in vcf_file.fetch():

	n_total=0
	n_hom_ref=0
	n_hom_alt=0
	n_het=0

	hom_ref=np.zeros(2, dtype=int)
	hom_alt=np.zeros(2, dtype=int)
	het=np.zeros(2, dtype=int)

	for sample, info in var.samples.items():
		if not all(info.alleles):
			continue


		a0=info.alleles[0]
		a1=info.alleles[1]

		if a0!=a1:
			n_het+=1
			het+=info["ARD"]
		else:
			if a0==var.ref:
				n_hom_ref+=1
				hom_ref+=info["ARD"]
			else:
				n_hom_alt+=1
				hom_alt+=info["ARD"]

		n_total+=1



	try:
		r=float(het[0])/float(np.sum(het))
	except ZeroDivisionError:
		r=np.nan


	if n_het<2 or np.sum(het)<50:
		p=np.nan
	else:
		p = binom_test(het, p=0.5)

	
	print '\t'.join(
		map(str,
			[
				var.contig, 
				var.start, 
				var.start+1, 
				var.id, 
				var.ref, 
				var.alts[0], 
				n_total, 
				n_hom_ref, 
				':'.join(map(str, hom_ref)), 
				n_hom_alt, 
				':'.join(map(str, hom_alt)), 
				n_het, 
				':'.join(map(str, het)),
				r,
				p
			])
		)
