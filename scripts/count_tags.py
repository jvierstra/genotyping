#!/bin/env python

import pysam
import sys
import logging

import numpy as np

logging.basicConfig(stream = sys.stdout, level = 20)


#sample = "/net/seq/data/aggregations/LN20963/aggregation-3871/LN20963.GRCh38_no_alts.sorted.bam"

vcf_file = sys.argv[1]
sample = sys.argv[2]
#bam_file= sys.argv[3]

vcf = pysam.VariantFile(vcf_file, "rb")
sam = pysam.AlignmentFile(sample, "rb" )

class GenotypeError(Exception):
	pass

class ReadError(Exception):
	
	ERROR_ALIGNMENT = (0, "Read alignment problematic (QC fail, duplicate, or MAPQ < 1)")
	ERROR_5PROXIMITY = (1, "Variant too close to 5' end of tag") 
	ERROR_BASEQ = (2, "Base quality < 20")
	ERROR_GENOTYPE = (3, "Base does not match reference or expected alternate allele")
	ERROR_MISMATCH = (4, "Read contains too many mismatches")
	
	def __init__(self, e):
		self.value = e[0]
		self.message = e[1]

#region="chr8:6648942-6648943"
#for rec in vcf.fetch(region = region):
for rec in vcf.fetch():

	#print rec.contig, rec.start, rec.ref, rec.ref

	ref = rec.ref
	alt = rec.alts[0]

	n_ref = 0
	n_alt = 0
	n_non = np.zeros(5, dtype = int)

	try:

		genotype = rec.samples[sample].alleles

		if not all(genotype):
			raise GenotypeError("No genotype called at site") # not genotyped

		for read in sam.fetch(rec.contig, rec.start, rec.start+1):
		  	
		  	try:

			  	if read.is_qcfail or read.is_duplicate or read.mapping_quality < 1:
			  		raise ReadError(ReadError.ERROR_ALIGNMENT)

				pos = read.reference_start
				offset = rec.start - pos

				nuc = read.query_sequence[offset]			
				qual = read.query_qualities[offset] #TODO: FILTER BY QUALITY AND BY # MISMATCHES

				mm = int(read.get_tag("XM", with_value_type = False))

				offset_5p = read.reference_end-1-rec.start if read.is_reverse else rec.start-read.reference_start
				
				if offset_5p <= 3:
					raise ReadError(ReadError.ERROR_5PROXIMITY)

				if qual < 20:
					raise ReadError(ReadError.ERROR_BASEQ)

				if nuc == ref:
					# Edit distance max 1 if read maps to reference
					if mm > 1:
						 raise ReadError(ReadError.ERROR_MISMATCH)
					n_ref += 1
				elif nuc == alt:
					if mm > 2:
						raise ReadError(ReadError.ERROR_MISMATCH)
					n_alt +=1 
				else:
					raise ReadError(ReadError.ERROR_GENOTYPE)

			except ReadError, e:
				
				n_non[e.value] += 1

				logging.debug(read.query_name + ": " + e.message)
				continue

			except IndexError, e:
				
				logging.debug(read.query_name + ": " + e.message)
				continue

#			except:
#				pass

	except GenotypeError, e:
		genotype = ('.', '.')
		
		logging.debug(e)
		pass

	except KeyboardInterrupt:
		sys.exit(1)

#	except:
#		pass
	
	print "%s\t%d\t%d\t%s\t%s:%d:%d\t%s" % (rec.contig, rec.start, rec.start+1, ref+"/"+alt, '/'.join(genotype), n_ref, n_alt, ':'.join(map(str, n_non)))
			
sam.close()
vcf.close()
