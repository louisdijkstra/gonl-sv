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
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import pyplot
from random import choice

import plfit
import plplot


# Global variables
LD_THRES = .8 	# R^2 threshold for a deletion and SNP to be considered 'related'
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



		
			

def sample():
	results = [] 

	chromosomes = range(1,23) 

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

		
		

	

def printSNPData (data):
	for snp in data:
		print snp


# Reads SNP-Deletion data
# 
# Args:
#	chr		- number of chromosome (1-22)
#	hit 		- 'all', 'hit', or 'nonhit'
#	type_snps 	- 'intergenic', 'intronic' or 'exonic' 
def readSNPData (chr, hit='all', type_snps=['intergenic', 'intronic', 'exonic']):
	
	chr = str(chr)		
	
	file_snps = open('SNP_max_r2_chr' + chr + '.txt', 'r')
	
	data = [] ; # read data is stored here

	# Read through all lines...
	for line in file_snps.readlines():
		array = [value for value in line.split()]
		
		hit_snp 	= int(array[1])
		snp_pos 	= int(array[2])
		type_snp 	= array[3]
		snp_maf 	= float(array[4])
		cnv_pos 	= int(array[5])
		cnv_length 	= int(array[6])
		r2_max 		= float(array[7])
		p 		= float(array[8])
		dist_snp_cnv 	= int(array[9])
		dist_tss 	= int(array[10])

		if hit == 'hit' and hit_snp == 1:
			if type_snp in type_snps:			
				data.append([hit_snp, snp_pos, type_snp, snp_maf, cnv_pos, cnv_length, r2_max, p, dist_snp_cnv, dist_tss])					
		if hit == 'nonhit' and hit_snp == 0:
			if type_snp in type_snps:			
				data.append([hit_snp, snp_pos, type_snp, snp_maf, cnv_pos, cnv_length, r2_max, p, dist_snp_cnv, dist_tss])					
		if hit == 'all':
			if type_snp in type_snps:		
				data.append([hit_snp, snp_pos, type_snp, snp_maf, cnv_pos, cnv_length, r2_max, p, dist_snp_cnv, dist_tss])					

	file_snps.close()			
	return data 	







def meanMaxR2SNPs (type_snps=['intergenic', 'intronic', 'exonic']):
	mean = 0 
	n = 0 

	for chr in range(1,23):
		data = readSNPData(chr, 'nonhit', type_snps)
		for snp in data:
			if not math.isnan(snp[6]):
				mean = mean + snp[6] 
				n = n + 1

	return mean / float(n)



def meanMaxR2HitSNPs (type_snps=['intergenic', 'intronic', 'exonic']):
	mean = 0 
	n = 0 

	for chr in range(1,23):
		data = readSNPData(chr, 'hit', type_snps)
		for snp in data:
			if not math.isnan(snp[6]):
				mean = mean + snp[6] 
				n = n + 1

	return mean / float(n)


def meanMaxR2HitSNPsPerChromosome (chr, type_snps=['intergenic, intronic, exonic']):
	data = readSNPData(chr,'hit',type_snps)
	mean = 0 
	for snp in data:
		if not math.isnan(snp[6]):
			mean = mean + snp[6]
		else:
			print 'r^2 is nan'

	return mean / len(data)


def overlapSNPs():
	for chr in range(1,23):
		n_hit = 0
		n_nonhit = 0  
		data = readSNPData(chr)
		for snp in data:
			if snp[8] == -1:
				if snp[0] == 1:
					n_hit = n_hit + 1
				else:
					n_nonhit = n_nonhit + 1
		print 'Chromosome ' + str(chr) + ':', n_hit, n_nonhit
	
def plotDistanceToTSS ():
	dtss_hit = [] 
	dtss_nonhit = [] 
	
	for chr in range(1,23):
		data = readSNPData(chr, hit='all')
		for snp in data:
			if snp[0] == 1: # hit-SNP
				dtss_hit.append(snp[9])
			else:
				dtss_nonhit.append(snp[9])
	
	

	dtss_hit.sort()
	dtss_nonhit.sort()
	pyplot.subplot(2,2,1)
	pyplot.plot(dtss_hit, color='red')
	#pyplot.xscale('log')
	#pyplot.yscale('log')
	pyplot.subplot(2,2,2)
	pyplot.plot(dtss_nonhit, color='blue')
	
	pyplot.subplot(2,2,3)
	pyplot.hist(dtss_hit,1000, facecolor='red', alpha=0.5)
	pyplot.xscale('log')
	pyplot.yscale('log')
	pyplot.subplot(2,2,4)
	pyplot.hist(dtss_nonhit,1000, facecolor='blue', alpha=0.5)
	pyplot.xscale('log')
	pyplot.yscale('log')
	pyplot.show() 

	#plt.hist(dtss_hit,1000, facecolor='red', alpha=0.5)
	#plt.hist(dtss_nonhit,1000, facecolor='red', alpha=0.5)
	#plt.xlabel('Distance to TSS')
	#plt.ylabel('Probability')
	#plt.title('Distance to TSS Distribution (GWAS/non-GWAS SNPs)')

	#plt.show() 


def plotMaxR2 ():
	r2s_hit = [] 
	r2s_nonhit = [] 
	for chr in range(1,23):
		data = readSNPData(chr, hit='hit')
		for hitsnp in data:
			if not math.isnan(hitsnp[6]):
				r2s_hit.append(hitsnp[6])
		data = readSNPData(chr, hit='nonhit')
		print 'chromosome ', chr
		for snp in data:
			if not math.isnan(snp[6]):
				r2s_nonhit.append(snp[6])
	
	bins = np.linspace(0.0,1.0,100)
	w = np.ones_like(r2s_hit) / len(r2s_hit)	
	plt.hist(r2s_hit, bins, weights=w, facecolor='red', alpha=0.5)
	w = np.ones_like(r2s_nonhit) / len(r2s_nonhit)
	plt.hist(r2s_nonhit, bins, weights=w, facecolor='blue', alpha=0.5)
	plt.xlabel(r'$R^2$')
	plt.ylabel('Probability')
	plt.title('SNP-Deletion Association')

	plt.show() 


#print meanMaxR2SNPs()
#print meanMaxR2SNPs('intergenic')
#print meanMaxR2SNPs('intronic')
#print meanMaxR2SNPs('exonic')

# overlapSNPs()

#plotMaxR2()

# sample()

plotDistanceToTSS()
