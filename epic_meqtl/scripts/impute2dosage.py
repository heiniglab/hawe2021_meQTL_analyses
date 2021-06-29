#!/usr/bin/env python

# read the impute output and transform to allele dosage
# first apply a filter to the genotype posteriors
# then apply a filter on the efficacy (number of people genotyped per SNP)
# finally filter by minor allele freq

import sys
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
#parser.add_argument('file', metavar='input_file', help='input file')
parser.add_argument('-p', metavar='min_posterior', type=float, default=0.95,
                   help='minimal genotype posterior probability')
parser.add_argument('-e', metavar='min_efficacy', type=float, default=0.95,
                   help='minimal imputation efficacy')
parser.add_argument('-m', metavar='min_MAF', type=float, default=0.05,
                   help='minimal minor allele frequency')
parser.add_argument('-o', metavar='output_file', default=None,
                   help='output file')
args = parser.parse_args()

cutoff = args.p
eff_cutoff = args.e
MAF_cutoff = args.m
outfile = sys.stdout
if (args.o):
    outfile = open(args.o, "w")

# file format is id1, id2, position, allele1, allele2, pers1 [P(gt1),P(gt2),]
#f = open(args.file)
# workaround to allow stdin for now
f = sys.stdin;
line = f.readline()
while (True):
    items = line.split(" ")
    # just print out the first five items
    out = ("\t").join(items[0:5])
    nsamples = (len(items) - 5) // 3
    # then go through the samples
    efficacy = 0
    MAF = 0
    for i in range(nsamples):
        # check if any of the three probabilities is above the treshold
        ok = False
        dosage = 0
        for gt in range(3):
            posterior = float(items[5 + gt + 3 * i])
            dosage = dosage + gt * posterior
            if (posterior > cutoff):
                efficacy = efficacy + 1
                ok = True
        MAF = MAF + dosage
        # if (not ok):
        #     dosage = "NA"
        out = "%s\t%s" % (out, dosage)
    efficacy = float(efficacy) / nsamples
    MAF = float(MAF) / (2 * nsamples)
    if (MAF > 0.5):
        MAF = 1 - MAF
#    print("MAF:" + str(MAF) + " ; Efficacy:" + str(efficacy))
#    if (efficacy > eff_cutoff and MAF > MAF_cutoff):
    if (MAF > MAF_cutoff):
        outfile.write(out + "\n")
    line = f.readline()
    if (not line):
        break


