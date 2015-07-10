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
from BenjaminiHochbergFDR import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] <snp-deletion-pair-file1> ... <snp-deletion-pair-file_m>

Reads in multiple '.pairs' files and returns the percentage of pairs deemed significant. 
The option '--list' returns the list of all significant pairs. 
"""

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("--list", action="store_true", dest="output_list", default=False, 
					help = "List of significant pairs is outputted.")
	parser.add_option("-q", action="store", dest="q", default=0.05, type=float,
				help = "Q-value for Benjaminin-Hochberg's FDR procedure. (Default = 0.05)") 
	parser.add_option("-r", action="store", dest="r2_threshold", default=0.8, type=float,
				help = "R^2 threshold for determining material significance. (Default = 0.8)") 
	(options, args) = parser.parse_args()

	if (len(args) < 1):
		parser.print_help()
		return 1

	r = []
	p_values = [] 

	n_pairs 		= 0  # total number of pairs
	n_mat_sign_pairs 	= 0  # materially significant pairs 
	n_stat_mat_sign_pairs 	= 0  # materially and statistically significant pairs
		
	for filename in args: 
		pairs_file = open(filename, 'r')
		
		for line in pairs_file:
			if line[0] != '-':
				continue
			values = line.split() 
			n_pairs += 1 
			r.append(float(values[3])) 
			p_values.append(float(values[4]))
			if r[-1] ** 2.0 >= options.r2_threshold:
				n_mat_sign_pairs += 1 
			
	r_stat_sign = benjamini_hochberg(options.q, p_values, r) ; 
	for r_value in r_stat_sign:
		if r_value ** 2 >= options.r2_threshold:
			n_stat_mat_sign_pairs += 1

	print("# SNP-deletion pairs: ", n_pairs)
	print("# materially significant SNP-deletion pairs (R^2 >= %f): "%options.r2_threshold, n_mat_sign_pairs)
	print("# statistically significant SNP-deletion pairs (q = %f): "%options.q, len(r_stat_sign))
	print("# materially/statistically significant SNP-deletion pairs: ", n_stat_mat_sign_pairs)

	if options.output_list: 
		pass
		


if __name__ == '__main__':
	sys.exit(main())

