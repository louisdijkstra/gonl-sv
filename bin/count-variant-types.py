#!/usr/bin/env python

"""
Copyright (C) 2015 Louis Dijkstra

This file is part of gonl-sv

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import math
import vcf

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))[:-3] + 'python')
from Variant import *

__author__ = "Louis Dijkstra"

usage = """%prog <haplotype-file1> <haplotype-file2> ... <haplotype-file_m>

Counts the number of variant types (deletion, insertion, 
GWAS SNP, non-GWAS SNP) in all given haplotype files. 
"""

def main():
	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()

	if (len(args) < 1):
		parser.print_help()
		return 1

	n_deletions 	= 0
	n_insertions 	= 0 
	n_gwas_snps 	= 0 
	n_nongwas_snps 	= 0 
	n_other	= 0 

	for filename in args:
		haplo_file = open(filename, 'r')
		for line in haplo_file:
			variant_type = line.split()[1][0]
			if variant_type == '*':
				n_nongwas_snps += 1 
			elif variant_type == '-':
				n_deletions += 1 
			elif variant_type == '+':
				n_insertions += 1 
			elif variant_type == '!':
				n_gwas_snps += 1 
			else:
				n_other += 1

	print("Number of SNPs: ", n_gwas_snps + n_nongwas_snps) 
	print("   # GWAS SNPs ('!'):\t\t", n_gwas_snps)
	print("   # non-GWAS SNPs ('*'):\t", n_nongwas_snps)	
	print("")
	print("Number of Indels: ", n_deletions + n_insertions) 
	print("   # deletions ('-'):\t\t", n_deletions)
	print("   # insertions ('+'):\t\t", n_insertions)	
	print("\nNumber of unspecified variants:\t", n_other) 
	print("\nTOTAL: ", n_gwas_snps + n_nongwas_snps + n_deletions + n_insertions + n_other)
if __name__ == '__main__':
	sys.exit(main())

