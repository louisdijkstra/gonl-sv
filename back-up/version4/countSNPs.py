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
import scipy as sp
import scipy.stats
from scipy.stats.stats import pearsonr

from SNPReader import *

__author__ = "Louis Dijkstra"

usage = """%prog <threshold>

Counts the number of GWAS SNPs and non-GWAS SNPs with a minor allele frequency above a given threshold 
"""	

def main():

	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()

	if (len(args)!=1):
		parser.print_help()
		return 1


	n_hit_snps = 0
	n_nonhit_snps = 0 
	for autosome in range(1,23):
		print("Reading autosome %d"%autosome)
		snp_reader = SNPReader(autosome, max_major_af = 1.0 - float(args[0]))
		snp_reader.read()
		
		for snp in snp_reader.snps:
			if snp.isHitSNP():
				n_hit_snps += 1
			else:
				n_nonhit_snps += 1
		print("# hit SNPs: %d\t# non-hit SNPs: %d"%(n_hit_snps, n_nonhit_snps))

if __name__ == '__main__':
	sys.exit(main())
