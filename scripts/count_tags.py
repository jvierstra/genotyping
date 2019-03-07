#!/bin/env python
# Jeff Vierstra 2018
# TODO:
# --add filters/etc. as option

import pysam
import sys
import logging
import argparse

from argparse import ArgumentParser


import numpy as np

logging.basicConfig(stream = sys.stdout, level = 20)

def parse_options(args):

    parser = argparse.ArgumentParser(description = "Count tags by allele")

    parser.add_argument("vcf_file", metavar = "vcf_file", type = str,
                        help = "Path to GZVCF-format tag sequence file")

    parser.add_argument("bam_file", metavar = "bam_file", type = str, 
                        help = "Path to BAM-format tag sequence file")

    parser.add_argument("sample", metavar = "sample", type = str, 
                        help = "Sample name (must be in the VFC file)")

    return parser.parse_args(args)

class GenotypeError(Exception):
    pass

class MatePairError(Exception):
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


def main(argv = sys.argv[1:]):

    args = parse_options(argv)

    vcf = pysam.VariantFile(args.vcf_file, "rb")
    sam = pysam.AlignmentFile(args.bam_file, "rb" )
    
    seen_mate_pairs = set()

    for rec in vcf.fetch():

        #print rec.contig, rec.start, rec.ref, rec.ref

        ref = rec.ref
        alt = rec.alts[0]

        n_ref = 0
        n_alt = 0
        n_non = np.zeros(5, dtype = int)

        seen_mate_pairs.clear()
        
        try:

            genotype = rec.samples[args.sample].alleles
            #ff = dict((x, y) for x, y in rec.samples[args.sample].items())

            # Check whether genotyping is of high quality and there is enough reads to support it

            if not all(genotype):
                raise GenotypeError("No genotype called at site for this individual") # not genotyped

            # Go into BAM file and get the reads

            for read in sam.fetch(rec.contig, rec.start, rec.start+1):
                try:

                    if read.query_name in seen_mate_pairs:
                        raise MatePairError("Mate pair already counted!")

                    if read.is_qcfail or read.is_duplicate or read.mapping_quality < 1:
                        raise ReadError(ReadError.ERROR_ALIGNMENT)

                    pos = read.reference_start
                    offset = rec.start - pos

                    nuc = read.query_sequence[offset]
                    qual = read.query_qualities[offset] #TODO: FILTER BY QUALITY AND BY # MISMATCHES

                    # XM tag stores number of mismatches
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
                        # Edit distance max 1 if read maps to alternate
                        if mm > 2:
                            raise ReadError(ReadError.ERROR_MISMATCH)
                        n_alt +=1 
                    else:
                        raise ReadError(ReadError.ERROR_GENOTYPE)

                    # a mate tag was successfully counted, exclude further mates from analysis
                    seen_mate_pairs.add(read.query_name)

                except ReadError, e:
                    
                    n_non[e.value] += 1

                    logging.debug(read.query_name + ": " + e.message)
                    continue

                except MatePairError, e:

                    logging.debug(read.query_name + ": " + e.message)
                    continue

                except IndexError, e:
                    
                    logging.debug(read.query_name + ": " + e.message)
                    continue
    
    #            except:
    #                pass

        except GenotypeError, e:
            genotype = ('.', '.')
            
            logging.debug(e)
            pass

        except KeyboardInterrupt:
            sys.exit(1)

    #    except:
    #        pass
        
        print "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%s" % (rec.contig, rec.start, rec.start+1, '/'.join(genotype), n_ref+n_alt, n_ref, n_alt, ':'.join(map(str, n_non)))
                
    sam.close()
    vcf.close()

    return 0
    
if __name__ == "__main__":
    sys.exit(main())

