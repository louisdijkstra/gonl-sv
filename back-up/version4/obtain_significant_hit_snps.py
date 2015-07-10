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

from SNPDeletionPair import *
from BenjaminiHochbergFDR import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] 

Outputs a list of hit SNP-deletion pairs that are deemed both statistically and materially significant. 
"""

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("-c", action="store", dest="chromosome", default=None, type=int,
				  		help="Only this autosome is taken into account. Must be an integer between 1 and 22.")
	parser.add_option("-q", action="store", dest="q", default=.05, type=float,
				  		help="Q-value used in the Hochberg & Benjamini FDR control procedure. (Default = 0.05)")
	parser.add_option("-r", action="store", dest="r2_threshold", default=.8, type=float,
				  		help="Minimal Pearson's R^2 level needed for a SNP-deletion pair to be considered materially significant. (Default = 0.8)")
	(options, args) = parser.parse_args()
	
	list_of_autosomes = range(1,23)
	if options.chromosome is not None:
		list_of_autosomes = [options.chromosome]

	# walk through all chromosomes
	hit_snp_deletion_pairs, p_values = [], [] # allocate memory
	for chromosome in list_of_autosomes: 
		results_file = open("Results/results.chr" + str(chromosome) + ".parents.txt")
		hit_pairs_on_chromosome = readSNPDeletionPairs(results_file, hit_snps_only = True, non_hit_snps_only = False)
		for pair in hit_pairs_on_chromosome:
			p_values.append(pair.p)
		hit_snp_deletion_pairs += hit_pairs_on_chromosome
		results_file.close()
				

	statistically_significant_pairs = benjamini_hochberg (options.q, p_values, hit_snp_deletion_pairs)	
	
	snp_deletion_pairs_of_interest = [] 
	for snp_deletion_pair in statistically_significant_pairs:
		if snp_deletion_pair.materiallySignificant(r2_threshold = options.r2_threshold):
			snp_deletion_pairs_of_interest.append(snp_deletion_pair)

	# output the results
	print('Initial number of hit SNP-deletion pairs considered: ', len(hit_snp_deletion_pairs))
	print('Number of hit SNP-deletion pairs considered statistically significant: ', len(statistically_significant_pairs))
	print('Number of hit SNP-deletion pairs deemed both statistically & materially significant: ', len(snp_deletion_pairs_of_interest))
	print('')
	print("CHR\tSNP_POS\tHITALLELE\tDEL_POS\tDEL_LENGTH\tR\tP\tA\tB\tC\tD")
	for snp_deletion_pair in snp_deletion_pairs_of_interest:
		snp_deletion_pair.print2()

if __name__ == '__main__':
	sys.exit(main())
