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

import math 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot
from random import choice


# Global variables
LD_THRES = .05 	# R^2 threshold for a deletion and SNP to be considered 'related'
DIST_TSS_THRES = 5000 
MAF_THRES = 0.10
N_SAMPLES = 100 # number of samples 



def sampleChr (chr): 
	# Read in data for chromosome
	data = readSNPData(chr)

	# walk through data - determine # of hit snps and 
	# discard any SNP that overlaps with its deletion
	n_hitsnps = 0 	
	lines_to_be_removed = [] 

	for i in range(len(data)):
		if data[i][8] == -1 or math.isnan(data[i][6]): # snp overlaps with deletion
			lines_to_be_removed.append(i)
		elif data[i][0] == 1:
			n_hitsnps = n_hitsnps + 1 

	for i in sorted(lines_to_be_removed, reverse=True):
		del data[i]

	
	control_indices = [] # contains indices of non-hit SNPs for each hit-SNP
	
	for hitsnp_index in range(n_hitsnps):
		similar_nonhitsnps = [] # list of indices of similar non-hit SNPs
		
		hit_type = data[hitsnp_index][2] # intergenic/intronic/exonic
		hit_maf = data[hitsnp_index][3]	 # minor allele frequency
		hit_dtss = data[hitsnp_index][9] # distance to TSS

		# runs over all non-hit snps and find the matching pairs
		for i in range(n_hitsnps, len(data)):
			# non-hit-SNP similar to hit-SNP?
			if hit_type == data[i][2] and abs(hit_maf - data[i][3]) <= MAF_THRES and abs(hit_dtss - data[i][9]) <= DIST_TSS_THRES:
				similar_nonhitsnps.append(i)

		# TODO: what if similar_nonhitsnps is empty? 
		if len(similar_nonhitsnps) == 0:
			print 'WARNING: no matching non-hit-SNP found. Constraints too stringent'
		
		#print len(similar_nonhitsnps)

		control_indices.append(similar_nonhitsnps)

	sample_results = [] # results for controls are stored here

	# start sampling 
	for sample in range(N_SAMPLES):
		
		control_LD = 0.0 # number of control SNPs in LD with a deletion
		for snp in range(n_hitsnps):
			control_snp = data[choice(control_indices[snp])]
			if control_snp[6] >= LD_THRES:
				control_LD = control_LD + 1 

		sample_results.append(control_LD / float(n_hitsnps))
		#print 'chr: ' + str(chr) + '| sample: ' + str(sample + 1)
	
	# return percentage of hit snps in LD,
	# average percentage of non-hit snps in LD, and
	# number of hit-snps in this chromosome
	mean_control = 0.0 
	for sample in range(N_SAMPLES):
		mean_control = mean_control + sample_results[sample]

	mean_control = mean_control / N_SAMPLES 

	case = 0.0 
	for hitsnp_index in range(n_hitsnps):
		if data[hitsnp_index][6] >= LD_THRES:
			case = case + 1 
		
	case = case / float(n_hitsnps)

	return [n_hitsnps, case, sample_results]	 



		
			

def sample(chromosomes):
	results = [] 

	for chr in chromosomes:
		r = sampleChr(chr)
		print str(chr) + '\t' + ' | # hit snps: ', r[0], '| % in LD: ', r[1] 
		results.append(r)

	perc_hit_LD = 0.0 
	n_total_hitsnps = 0.0 
	samples = [0.0] * N_SAMPLES 

	for i in range(len(chromosomes)):
		n_hitsnps = results[i][0]
		n_total_hitsnps = n_total_hitsnps + n_hitsnps
		perc_hit_LD = perc_hit_LD + n_hitsnps * results[i][1]
		for s in range(N_SAMPLES):
			samples[s] = samples[s] + n_hitsnps * results[i][2][s]

	perc_hit_LD = perc_hit_LD / n_total_hitsnps 
	for s in range(N_SAMPLES):
		samples[s] = samples[s] / n_total_hitsnps

	print perc_hit_LD
	print samples


def retrieveHitSNPPositions (chr):
	chr = str(chr)		
	file_snps = open('/ufs/dijkstra/Projects/SNPs_LD_deletions/SNP_max_r2_chr' + chr + '.txt', 'r')
	
	hitSNPpositions = [] 
	
	for line in file_snps.readlines():
		l = [value for value in line.split()]
		if int(l[1]) == 1:
			hitSNPpositions.append(int(l[2]))
		else:
			break 

	return hitSNPpositions 

def printSNPData(data):
	for snp in data:
		print snp 


# Reads SNP-Deletion data
# 
# Args:
#	chr		- number of chromosome (1-22)
#	type_snps 	- 'intergenic', 'intronic' or 'exonic' 
def readSNPData (chr, type_snps=['intergenic', 'intronic', 'exonic']):
	
	hitsnps = retrieveHitSNPPositions(chr)

	chr = str(chr)		
	file_snps = open('del_snp_pearson_r_chr' + chr + '.phased', 'r')
	
	data = [] ; # read data is stored here

	# Read through all lines...
	for line in file_snps.readlines():
		array = [value for value in line.split()]
			
		snp_pos 	= int(array[1])
		hitsnp 		= 0 
		if snp_pos in hitsnps:
			hitsnp = 1 
		
		type_snp 	= array[2]
		snp_maf 	= float(array[3])
		cnv_pos 	= int(array[4])
		cnv_length 	= int(array[5])
		cnv_maf 	= float(array[6])
		r2_max 		= float(array[7])
		p 		= float(array[8])
		dist_snp_cnv 	= int(array[9])
		dist_tss 	= int(array[10])

		if type_snp in type_snps:
			data.append([hitsnp, snp_pos, type_snp, snp_maf, cnv_pos, cnv_length, r2_max, p, dist_snp_cnv, dist_tss])

		if hitsnp == 1: # TODO remove
			print data[-1]

	file_snps.close()			
	return data 	

sample([1])



