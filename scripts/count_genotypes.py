import sys
import logging
import argparse

from argparse import ArgumentParser

import numpy as np
import scipy
from scipy.stats import binom_test

logging.basicConfig(stream = sys.stdout, level = 20)

def parse_options(args):

    parser = ArgumentParser(description = "Aggregate summary of genotype read depth and allelic ratios")

    parser.add_argument("--min_het_reads", metavar = "min_het_reads", type = int, default=50, 
                        help = "Mininum number of HETEROZYGOUS reads to test for imbalance")

    parser.add_argument("--min_het_samples", metavar = "min_het_samples", type = int, default=1, 
                        help = "Mininum number of HETEROZYGOUS samples to test for imbalance")

    parser.add_argument("--p", metavar = "p", type = float, default=0.5, 
                        help = "Binomial paramter p used in imbalance test (default: %(default)s)")

    return parser.parse_args(args)


def format_counts(a):
    return ':'.join(map(str, a))

def main(argv = sys.argv[1:]):

    args = parse_options(argv)

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

        if (heterozygous[0]+heterozygous[1] >= args.min_het_reads) and n_heterozygous>args.min_het_samples:
            p = binom_test(heterozygous[:2], p = args.p)
        else:
            p = np.nan

        print "\t".join(map(str, [chrom, start, end, "%s/%s" % (ref, alt) , total, n_homozygous_ref, format_counts(homozygous_ref), n_homozygous_alt, format_counts(homozygous_alt), n_heterozygous, format_counts(heterozygous), imbalance, p]))

if __name__ == "__main__":
    sys.exit(main())

