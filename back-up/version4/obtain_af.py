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
import random
import scipy as sp
import scipy.stats
from scipy.stats.stats import pearsonr

from DeletionReader import *
from SNPReader import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] 

Obtains and returns the allele frequency of all HumanOmni SNPs and deletions. The allele frequency 
for the SNP is with respect to the reference allele. The AF for the deletions is equal to the percentage 
of haplotypes carrying the indel. Results are based on the phase of the parents.
"""


def main():
	parser = OptionParser(usage=usage)
	parser.add_option("--deletions-only", action="store_true", dest="deletions_only", default=False,
						help="The allele frequency of deletions only.")
	parser.add_option("--snps-only", action="store_true", dest="snps_only", default=False, 
						help="The allele frequency of SNPs only.")
	parser.add_option("-c", action="store", dest="chromosome", default=None, type=int,
				  		help="Only this autosome is taken into account. Must be an integer between 1 and 22.")
	(options, args) = parser.parse_args()

	list_of_autosomes = range(1,23)
	if options.chromosome is not None:
		list_of_autosomes = [options.chromosome]

	
	deletion_reader = DeletionReader(max_major_af = 1.0) # used for reading in the deletions

	for autosome in list_of_autosomes:
		if not options.snps_only: # process the deletions as well: 
			deletion_reader.read(autosome)
			for deletion in deletion_reader.deletions:
				af = float(sum(deletion.haplotypes)) / float(len(deletion.haplotypes)) # determine allele frequency
				print(af)

		if not options.deletions_only: # process the snps as well: 
			snp_reader = SNPReader(autosome, max_major_af = 1.0)
			for snp in snp_reader.snps:
				af = float(sum(snp.haplotypes)) / float(len(snp.haplotypes))
				print(af)
		

if __name__ == '__main__':
	sys.exit(main())

