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

usage = """%prog <gwas-catalogue-file>

Reads in a sorts a GWAS catalogue file 
"""

def main():
	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()

	if (len(args)!=1):
		parser.print_help()
		return 1

	gwas_catalogue_file = open(args[0], 'r')

	snp_positions 	= [[] for i in range(1,23)] # one for every autosome
	hit_alleles 	= [[] for i in range(1,23)] # one for every autosome

	gwas_catalogue_file.next() # skip header

	for line in gwas_catalogue_file:
		values = line.split('\t')
		if len(values) != 4:
			continue
		autosome = returnAutosome(values[0].strip())
		if autosome < 23: 
			snp_positions[autosome - 1].append(int(values[1]))
			hit_alleles[autosome - 1].append(values[3]) 
			

	# sort array of positions per autosome
	for autosome in range(22):
		snp_positions[autosome], hit_alleles[autosome] = (list(t) for t in zip(*sorted(zip(snp_positions[autosome], hit_alleles[autosome]))))

	
	# print results to screen
	for autosome in range(22):
		for i in range(len(snp_positions[autosome])):
			print("%d %d %s"%(autosome+1, snp_positions[autosome][i], hit_alleles[autosome][i]), end = '')

if __name__ == '__main__':
	sys.exit(main())

