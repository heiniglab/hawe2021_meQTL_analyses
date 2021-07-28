#!/bin/python3.6

import sys
import argparse

parser = argparse.ArgumentParser(description='Adjust p-val sorted matrix QTL output for MT.')
parser.add_argument('-i', metavar='input_file', help='input file')
#parser.add_argument('-o', metavar='output_file', help='output file')
args = parser.parse_args()

#output = open(args.o, "w")

total_tests = 5738181750

with sys.stdout as output:
  last = None
  index = total_tests
  with open(args.i, "r") as inp:
    for line in inp:
      if(line.startswith("SNP")):
        continue;

      vals = line.split("\t")
      pv = float(vals[4])
      if(index == total_tests):
        bh = min(1, pv * (total_tests / index))
      else:
        bh = min(1, pv * (total_tests / index), last)

      output.write(line.strip() + "\t" + str(bh) + "\n")
      last = bh
      index -= 1
