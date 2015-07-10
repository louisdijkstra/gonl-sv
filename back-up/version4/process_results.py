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


from SNPDeletionPair import *
from BenjaminiHochbergFDR import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] 

Processes the SNP-Deletion pairs generated previously. It outputs the percentage of hit SNP-deletion pairs considered statistically significant. 
In addition, it outputs the percentage of matched non-hit SNP-deletion pairs considered statistically significant. 
"""

def determinePercentageSignificant(q, snp_deletion_pairs):
	"""Returns the percentage of SNP deletion pairs deemed significant after applying Benjamini Hochberg with the given q-value."""
	p_values = []
	for snp_deletion_pair in snp_deletion_pairs:
		p_values.append(snp_deletion_pair.p)
	significant_pairs = benjamini_hochberg(q, p_values, snp_deletion_pairs)
	return float(len(significant_pairs)) / float(len(snp_deletion_pairs))

def returnIndicesOfSimilarPairs(hit_snp_deletion_pairs, non_hit_snp_deletion_pairs, major_af_threshold = 0.05, distance_snp_del_threshold = 10000, distance_tss_threshold = 10000):
	all_indices = [] 
	for hit_snp_deletion_pair in hit_snp_deletion_pairs:
		indices = [] 
		for i in range(len(non_hit_snp_deletion_pairs)):
			if hit_snp_deletion_pair.similarTo(non_hit_snp_deletion_pairs[i], major_af_threshold = major_af_threshold, distance_snp_del_threshold = distance_snp_del_threshold, distance_tss_threshold = distance_tss_threshold):
				indices.append(i)
		all_indices.append(indices)
	return all_indices

def printList(l):
	if len(l) == 0:
		print('[]')
	else:
		print('[', end = '')
		for i in range(len(l) - 1):
			print(l[i], ', ', end = '')
		print(l[-1], ']')

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("-a", action="store", dest="major_af_threshold", default=.05, type=float,
				  		help="Maximal difference between major allele frequencies SNPs can have in order to be considered 'similar'. (Default = 0.05)")
	parser.add_option("-c", action="store", dest="chromosome", default=None, type=int,
				  		help="Only this autosome is taken into account. Must be an integer between 1 and 22.")
	parser.add_option("-d", action="store", dest="distance_snp_del_threshold", default=10000, type=int,
				  		help="Maximal distance between SNP and deletion for the pair to be considered. (Default = 10,000)")
	parser.add_option("-q", action="store", dest="q", default=.05, type=float,
				  		help="Q-value used in the Hochberg & Benjamini FDR control procedure. (Default = 0.05)")
	parser.add_option("-M", action="store", dest="number_of_samples", default=1000, type=int,
				  		help="Number of samples used. (Default = 1000)")
	parser.add_option("-t", action="store", dest="distance_tss_threshold", default=10000, type=int,
				  		help="Maximal difference in distance to the nearest TSS two SNP-deletion pairs might have in order to be considered 'similar'. (Default = 10,000)")
	(options, args) = parser.parse_args()

	list_of_autosomes = range(1,23) # autsomes that need to be processed
	if options.chromosome is not None: # option applies
		list_of_autosomes = [options.chromosome]

	hit_snp_deletion_pairs, non_hit_snp_deletion_pairs = [], [] # allocate memory
	indices_similar_pairs = [] # for every hit SNP a list of indices of non_hit snps that are similar
	n_non_hit_pairs = 0 
	# walk through all chromosomes
	for chromosome in list_of_autosomes:
		print('Preprocessing autosome', chromosome)
		results_file 		= open("Results/results.chr" + str(chromosome) + ".parents.txt", 'r')
		print('Reading in the SNP-deletion pairs from results.chr' + str(chromosome) + '.parents.txt')
		hit_snp_deletion_pairs_chr , non_hit_snp_deletion_pairs_chr = readHitAndNonHitSNPDeletionPairs(results_file)
		print('Done reading in the SNP-deletion pairs...\n')
		results_file.close()

		print('# hit-SNP deletion pairs\t:', len(hit_snp_deletion_pairs_chr))
		print('# non-hit-SNP deletion pairs\t:', len(non_hit_snp_deletion_pairs_chr), '\n')
		
		print('Find for every hit SNP-deletion pair every non-hit SNP-deletion pair that is similar.\n')
		indices_similar_pairs_chr 	= returnIndicesOfSimilarPairs(hit_snp_deletion_pairs_chr, non_hit_snp_deletion_pairs_chr, major_af_threshold = options.major_af_threshold, distance_snp_del_threshold = options.distance_snp_del_threshold, distance_tss_threshold = options.distance_tss_threshold)
		
		for x in indices_similar_pairs_chr:
			for i in range(len(x)):
				x[i] += n_non_hit_pairs
			
		indices_similar_pairs 		+= indices_similar_pairs_chr
		n_non_hit_pairs 		+= len(non_hit_snp_deletion_pairs_chr)
		hit_snp_deletion_pairs 		+= hit_snp_deletion_pairs_chr
		non_hit_snp_deletion_pairs	+= non_hit_snp_deletion_pairs_chr



	# Determine the percentage for the hit SNP deletion pairs
	print(determinePercentageSignificant(options.q, hit_snp_deletion_pairs))
	
	# Determine number of hit SNP deletion pairs that have not a single similar non-hit SNP equivalent
	no_similar_present = 0
	for i in range(len(hit_snp_deletion_pairs)):
		if len(indices_similar_pairs[i]) == 0:
			no_similar_present += 1

	print(no_similar_present, len(hit_snp_deletion_pairs))
	

	# Start with sampling!
	sampling_results = []
	for k in range(options.number_of_samples):
		sample = [None for i in range(len(hit_snp_deletion_pairs))] # sample of non-hit SNP-deletion pairs
		for i in range(len(hit_snp_deletion_pairs)):
			if len(indices_similar_pairs[i]) == 0: # if there are no similar pairs, select one at random
				sample[i] = random.choice(non_hit_snp_deletion_pairs) 
			else:
				sample[i] = non_hit_snp_deletion_pairs[random.choice(indices_similar_pairs[i])]
		sampling_results.append(determinePercentageSignificant(options.q, sample))
		print('sample ', k+1, ': ', sampling_results[-1])

	print(sampling_results)
			
		

if __name__ == '__main__':
	sys.exit(main())
