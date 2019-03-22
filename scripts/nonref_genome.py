#!/bin/env python

import sys
import shutil


import pyfaidx
import pysam



iupac = "XACMGRSVTWYHKDBN"
def get_iupac(x, y):
	i = iupac.find(x)
	j = iupac.find(y)
	return iupac[(i|j)]


orig_fasta_file=sys.argv[1]
vcf_file=sys.argv[2]
out_fasta_file=sys.argv[3]

print "copying original fasta file..."

shutil.copy2(orig_fasta_file, out_fasta_file)
shutil.copy2(orig_fasta_file+".fai", out_fasta_file+".fai")

print "making variant genome..."

fasta = pyfaidx.Fasta(out_fasta_file, mutable=True)
vcf = pysam.VariantFile(vcf_file)

for rec in vcf.fetch():
	
	if rec.filter.get("PASS") is None:	
		continue

	ref = rec.ref
	alt = rec.alts[0]

	if len(ref)>1 or len(alt)>1:
		continue

	ambig = get_iupac(ref, alt)
	fasta[rec.contig][rec.start] = ambig

	#print fasta[rec.contig][rec.start], ref, alt, ambig
