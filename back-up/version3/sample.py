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
import operator
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import pyplot
from random import choice
import scipy as sp
import scipy.stats

# Global variables
LD_THRES = .8 	# R^2 threshold for a deletion and SNP to be considered 'related'
D_TSS_THRES = 1000 
D_SNP_DEL_THRES = 100000
AF_THRES = 0.05
N_SAMPLES = 1000 # number of samples 



def snpsSimilar (hit_snp, nonhit_snp, af_thres = AF_THRES, d_tss_thres = D_TSS_THRES):
	if hit_snp['type'] == nonhit_snp['type'] and abs(hit_snp['d_tss'] - nonhit_snp['d_tss']) <= d_tss_thres:
		if abs(hit_snp['af'] - nonhit_snp['af']) <= af_thres or abs(hit_snp['af'] - (1 - nonhit_snp['af'])) <= af_thres:
			return True
	return False 


def returnControlIndices (hit_snps, nonhit_snps, af_thres = AF_THRES, d_tss_thres = D_TSS_THRES):
	control_indices = [] # will contain lists of indices. Each list contains indices for potential control snps
	
	for hit_snp in hit_snps:
		similar_nonhitsnps = [] 
		for i in range(len(nonhit_snps)):
			if snpsSimilar(hit_snp, nonhit_snps[i], af_thres, d_tss_thres):
				similar_nonhitsnps.append(i)
		
		if len(similar_nonhitsnps) == 0:
			print 'WARNING: no matching non-hit-SNP found. Constraints too stringent'
		
		control_indices.append(similar_nonhitsnps)

	return control_indices 

def snpsSimilar3 (hit_snp, nonhit_snp, af_thres = AF_THRES, d_tss_thres = D_TSS_THRES):
	if hit_snp['type'] == nonhit_snp['type'] and abs(hit_snp['d_tss'] - nonhit_snp['d_tss']) <= d_tss_thres:
		if abs(hit_snp['af'] - nonhit_snp['af']) <= af_thres or abs(hit_snp['af'] - (1 - nonhit_snp['af'])) <= af_thres:
			if abs(hit_snp['d_snp_del'] - nonhit_snp['d_snp_del']) <= D_SNP_DEL_THRES:
				return True
	return False 


def returnControlIndices3 (hit_snps, nonhit_snps, af_thres = AF_THRES, d_tss_thres = D_TSS_THRES):
	control_indices = [] # will contain lists of indices. Each list contains indices for potential control snps
	
	for hit_snp in hit_snps:
		similar_nonhitsnps = [] 
		for i in range(len(nonhit_snps)):
			if snpsSimilar3(hit_snp, nonhit_snps[i], af_thres, d_tss_thres):
				similar_nonhitsnps.append(i)
		
		if len(similar_nonhitsnps) == 0:
			print 'WARNING: no matching non-hit-SNP found. Constraints too stringent'
		
		control_indices.append(similar_nonhitsnps)

	return control_indices 

def sampleChromosome3 (chromosome, ld_thres = LD_THRES, d_tss_thres = D_TSS_THRES, af_thres = AF_THRES, n_samples = N_SAMPLES):
	[hit_snps, nonhit_snps] = readSNPData(chromosome)
	n_hit_snps = len(hit_snps)

	# control_indices contain lists of indices. Each list contains indices for potential control snps
	control_indices = returnControlIndices (hit_snps, nonhit_snps, af_thres, d_tss_thres)
	
	sample_results = [0.0] * n_samples 
	for sample in range(n_samples):
		for i in range(n_hit_snps):
			if len(control_indices[i]) != 0:
				hit_snp = hit_snps[i]				
				nonhit_snp = nonhit_snps[choice(control_indices[i])] # randomly select a similar nonhit SNP
				if abs(hit_snp['af'] - nonhit_snp['af']) <= abs(hit_snp['af'] - (1 - nonhit_snp['af'])):
					if nonhit_snp['r'] >= .84:
						sample_results[sample] = sample_results[sample] + 1
				else:
					if -1 * nonhit_snp['r'] >= .84:
						sample_results[sample] = sample_results[sample] + 1

	
	case = 0.0 
	for hit_snp in hit_snps:
		if hit_snp['r'] >= .84: 
			case = case + 1

	return [n_hit_snps, case, sample_results]


def plotSampleResults3 (case, control, n_samples = N_SAMPLES):
	bins = np.linspace(min(control),max(control),50)
	w = np.ones_like(control) / len(control)
	plt.hist(control, bins, weights=w, facecolor='blue', alpha=1)
	plt.axvline(x=case, color='r', linewidth=4)
	plt.xlabel(r'Percentage of SNPs with $R^2 > \ 0.7$')
	plt.ylabel('Probability')
	plt.show() 

def sample3(ld_thres = LD_THRES, d_tss_thres = D_TSS_THRES, af_thres = AF_THRES, n_samples = N_SAMPLES):
	sample_results = [0.0] * n_samples
	case = 0.0 
	n_hit_snps = 0.0 

	for chromosome in range(1,23):
		print 'Processing chromosome', chromosome
		[hits, case_chr, sample_chr] = sampleChromosome3(chromosome, ld_thres, d_tss_thres, af_thres, n_samples)
		case = case + case_chr
		n_hit_snps = n_hit_snps + hits 
		sample_results = map(operator.add, sample_results, sample_chr)		
		print 'Chr ' + str(chromosome) + '\t', hits, case_chr / hits
	
	for i in range(n_samples):
		sample_results[i] = sample_results[i] / n_hit_snps 
	
	print '\nOVERALL: ', n_hit_snps, case / n_hit_snps, sample_results  
	plotSampleResults3 (case / n_hit_snps, sample_results)


def sampleChromosome2 (chromosome, ld_thres = LD_THRES, d_tss_thres = D_TSS_THRES, af_thres = AF_THRES, n_samples = N_SAMPLES):
	[hit_snps, nonhit_snps] = readSNPData(chromosome)
	n_hit_snps = len(hit_snps)

	# control_indices contain lists of indices. Each list contains indices for potential control snps
	control_indices = returnControlIndices (hit_snps, nonhit_snps, af_thres, d_tss_thres)
	
	sample_results = [0.0] * n_samples 
	for sample in range(n_samples):
		for i in range(n_hit_snps):
			if len(control_indices[i]) != 0:
				hit_snp = hit_snps[i]				
				nonhit_snp = nonhit_snps[choice(control_indices[i])] # randomly select a similar nonhit SNP
				if abs(hit_snp['af'] - nonhit_snp['af']) <= abs(hit_snp['af'] - (1 - nonhit_snp['af'])):
					if nonhit_snp['r'] >= 0:
						sample_results[sample] = sample_results[sample] + 1
				else:
					if -1 * nonhit_snp['r'] >= 0:
						sample_results[sample] = sample_results[sample] + 1

	
	case = 0.0 
	for hit_snp in hit_snps:
		if hit_snp['r'] >= 0: 
			case = case + 1

	return [n_hit_snps, case, sample_results]


def sample2(ld_thres = LD_THRES, d_tss_thres = D_TSS_THRES, af_thres = AF_THRES, n_samples = N_SAMPLES):
	sample_results = [0.0] * n_samples
	case = 0.0 
	n_hit_snps = 0.0 

	for chromosome in range(1,23):
		print 'Processing chromosome', chromosome
		[hits, case_chr, sample_chr] = sampleChromosome2(chromosome, ld_thres, d_tss_thres, af_thres, n_samples)
		case = case + case_chr
		n_hit_snps = n_hit_snps + hits 
		sample_results = map(operator.add, sample_results, sample_chr)		
		print 'Chr ' + str(chromosome) + '\t', hits, case_chr / hits
	
	for i in range(n_samples):
		sample_results[i] = sample_results[i] / n_hit_snps 
	
	print '\nOVERALL: ', n_hit_snps, case / n_hit_snps, sample_results  
	plotSampleResults (case / n_hit_snps, sample_results)
	



def sampleChromosome (chromosome, ld_thres = LD_THRES, d_tss_thres = D_TSS_THRES, af_thres = AF_THRES, n_samples = N_SAMPLES):
	[hit_snps, nonhit_snps] = readSNPData(chromosome)
	n_hit_snps = len(hit_snps)

	# control_indices contain lists of indices. Each list contains indices for potential control snps
	control_indices = returnControlIndices (hit_snps, nonhit_snps, af_thres, d_tss_thres)
	
	sample_results = [0.0] * n_samples 
	for sample in range(n_samples):
		for i in range(n_hit_snps):
			if len(control_indices[i]) != 0:
				nonhit_snp = nonhit_snps[choice(control_indices[i])] # randomly select a similar nonhit SNP
				if nonhit_snp['r']**2 >= ld_thres:
					sample_results[sample] = sample_results[sample] + 1
	
	case = 0.0 
	for hit_snp in hit_snps:
		if hit_snp['r']**2 >= ld_thres: 
			case = case + 1

	return [n_hit_snps, case, sample_results]


def plotSampleResults (case, control, n_samples = N_SAMPLES):
	bins = np.linspace(min(control),max(control),50)
	w = np.ones_like(control) / len(control)
	plt.hist(control, bins, weights=w, facecolor='blue', alpha=1)
	plt.axvline(x=case, color='r', linewidth=4)
	plt.xlabel(r'Percentage of SNPs with $R > 0$')
	plt.ylabel('Probability')
	plt.show() 
	

def sample(ld_thres = LD_THRES, d_tss_thres = D_TSS_THRES, af_thres = AF_THRES, n_samples = N_SAMPLES):
	sample_results = [0.0] * n_samples
	case = 0.0 
	n_hit_snps = 0.0 

	for chromosome in range(1,23):
		print 'Processing chromosome', chromosome
		[hits, case_chr, sample_chr] = sampleChromosome(chromosome, ld_thres, d_tss_thres, af_thres, n_samples)
		case = case + case_chr
		n_hit_snps = n_hit_snps + hits 
		sample_results = map(operator.add, sample_results, sample_chr)		
		print 'Chr ' + str(chromosome) + '\t', hits, case_chr / hits
	
	for i in range(n_samples):
		sample_results[i] = sample_results[i] / n_hit_snps
	
	print '\nOVERALL: ', n_hit_snps, case / n_hit_snps, sample_results  
	plotSampleResults (case / n_hit_snps, sample_results)
	


def returnHitSNPPositions (chromosome):
	file_snps = open('/ufs/dijkstra/Projects/SNPs_LD_deletions/SNP_max_r2_chr' + str(chromosome) + '.txt', 'r')
	positions = [] # hit-SNP positions are stored here
	for snp in file_snps.readlines():
		data = [value for value in snp.split()]
		if int(data[1]) == 1:
			positions.append(int(data[2]))
		else: 
			break 
	file_snps.close()
	return positions 
















# Writes arrays to a prescripted output file 
def writeToOutputFile(output, chr, individuals):
	filename = None
	if 'children' in individuals:
		filename = 'chr' + str(chr) + '_children_metrics.txt' 
	else:
		filename = 'chr' + str(chr) + '_parents_metrics.txt'  
	output_file = open(filename, 'a')
	for line in output:
		for item in line:
			output_file.write(str(item) + '\t')
		output_file.write('\n')

	output_file.close() 



# Returns the distance between the SNP and the deletion (CNV)
# Returns 0 when the SNP lies within the deletion
def distanceSNP_DEL (snp_pos, del_start, del_length):
	if snp_pos < del_start:
		return del_start - snp_pos
	elif snp_pos > del_start + del_length:
		return snp_pos - del_start + del_length 
	else: # snp overlaps with the deletion
		return 0 

def readSNPDeletionPairs(chromosome, individuals):
	hitsnp_del_data = [] # output is stored here
	nonhitsnp_del_data = []  
	file_snps = None

	if 'children' in individuals:
		file_snps = open('chr' + str(chromosome) + '_children.txt', 'r')
	else: 
		file_snps = open('chr' + str(chromosome) + '_children.txt', 'r')

	for line in file_snps.readlines():
		l = [value for value in line.split()]
		snp_del = dict([('pos', int(l[1])), ('type', int(l[2])), ('dist_tss', int(l[3])), ('del_pos', int(l[4])), ('del_length', int(l[5])), ('a', int(l[6])), ('b', int(l[7])), ('c', int(l[8]))])
		if l[0] == '0':
			nonhitsnp_del_data.append(snp_del)
		else:
			hitsnp_del_data.append(snp_del)
	
	return [hitsnp_del_data, nonhitsnp_del_data]
	

def determineGamma_N(a,b,c,N):
	a = float(a)
	b = float(b)
	c = float(c)
	gamma = a / (a + c) - (a + b) / N 
	if gamma <= 0:
		return gamma * N / (a + b)
	elif a+c <= a + b:
		return gamma / (1 - (a+b)/N)
	else:
		return gamma / ((a + b)/(a+c) - (a+b)/N)
	

def determineR(a,b,c,N):
	d = N - a - b -c 
	return float((a * d - b * c)) / math.sqrt((a+b)*(a+c)*(b+d)*(c+d))

def fisher(a,b,c,N):
	table = [[a,c], [b, N - a - b - c]]

	[odds_ratio, p2] = sp.stats.fisher_exact(table, alternative='two-sided')
	[odds_ratio, pg] = sp.stats.fisher_exact(table, alternative='greater') 
	[odds_ratio, pl] = sp.stats.fisher_exact(table, alternative='less') 
	return [p2,pg,pl]

def compute(chromosome, individuals):
	N = 470
	if 'parents' in individuals:
		N = None # TODO!

	[hitsnp_del, nonhitsnp_del] = readSNPDeletionPairs(chromosome, individuals)
	output = [] 	
	for snp in hitsnp_del:
		dist_snp_del = distanceSNP_DEL(snp['pos'], snp['del_pos'], snp['del_length'])
		r = determineR(snp['a'], snp['b'], snp['c'], N)
		gamma_n = determineGamma_N(snp['a'], snp['b'], snp['c'], N)
		[p2,pg,pl] = fisher(snp['a'], snp['b'], snp['c'], N)
		output.append([1, snp['type'], snp['dist_tss'], dist_snp_del, r, gamma_n, p2,pg,pl])
		#print [1, snp['type'], snp['dist_tss'], dist_snp_del, r, gamma_n, p2,pg,pl]	
	for snp in nonhitsnp_del:
		dist_snp_del = distanceSNP_DEL(snp['pos'], snp['del_pos'], snp['del_length'])
		r = determineR(snp['a'], snp['b'], snp['c'], N)
		gamma_n = determineGamma_N(snp['a'], snp['b'], snp['c'], N)
		[p2,pg, pl] = fisher(snp['a'], snp['b'], snp['c'], N)
		output.append([0, snp['type'], snp['dist_tss'], dist_snp_del, r, gamma_n, p2,pg,pl])
		#print [0, snp['type'], snp['dist_tss'], dist_snp_del, r, gamma_n, p2,pg,pl]	
	writeToOutputFile(output, chromosome, individuals)
	

def readSNPData (chromosome):
	hit_snp_positions = returnHitSNPPositions (chromosome)
	file_snps = open('chr' + str(chromosome) + '.txt', 'r')
	hit_snps = [] # hit-snp data is stored here
	nonhit_snps = [] # non hit-snp data is stored here
	# Read through all lines
	for line in file_snps.readlines():
		l = [value for value in line.split()]		
		snp = [('hit', int(l[1])),('pos', int(l[2])), ('type', l[3]), ('af', float(l[4]))]
		snp.append(('del_pos', int(l[5])))
		snp.append(('del_length', int(l[6])))
		snp.append(('del_af', float(l[7])))
		snp.append(('r', float(l[8])))
		snp.append(('p', float(l[9])))
		snp.append(('d_snp_del', int(l[10])))
		snp.append(('d_tss', int(l[11])))
		snp.append(('min_dist_snp_del', int(l[12])))
		snp.append(('g', float(l[13])))
		snp.append(('g_n', float(l[14])))
		snp.append(('a', int(l[15]))) 
		snp.append(('b', int(l[16])))
		snp.append(('c', int(l[17])))
		snp.append(('d', int(l[18])))

		if int(l[2]) in hit_snp_positions:
			hit_snps.append(dict(snp))
		else:
			nonhit_snps.append(dict(snp))

	file_snps.close()
	return [hit_snps, nonhit_snps]



def plotRDistribution_nonhit ():
	nonhit_r = [] 
	for chromosome in range(1,23):
		print 'Processing chromosome', chromosome
		[hit_snps, nonhit_snps] = readSNPData(chromosome)
		for snp in nonhit_snps:
			nonhit_r.append(snp['r'])
			#print snp['r']
	bins = np.linspace(-1.0,1.0,100)
	w = np.ones_like(nonhit_r) / len(nonhit_r)
	plt.hist(nonhit_r, bins, weights=w, facecolor='blue', alpha=1)
	plt.xlabel(r'$R$')
	plt.ylabel('Probability')
	plt.title(r'non-GWAS SNP - Deletion combinations')
	plt.savefig('R_nonhit.png')
	plt.show() 
	


def readSNPData_Androniki ():
	hit_snp_positions = returnHitSNPPositions (20)
	file_snps = open('androniki_r_chr20.phased', 'r')
	hit_snps = [] # hit-snp data is stored here
	nonhit_snps = [] # non hit-snp data is stored here
	# Read through all lines
	for line in file_snps.readlines():
		l = [value for value in line.split()]		
		snp = [('hit', int(l[1])),('pos', int(l[2])), ('type', l[3]), ('af', float(l[4]))]
		snp.append(('del_pos', int(l[5])))
		snp.append(('del_length', int(l[6])))
		snp.append(('del_af', float(l[7])))
		snp.append(('r', float(l[8])))
		snp.append(('p', float(l[9])))
		snp.append(('d_snp_del', int(l[10])))
		snp.append(('d_tss', int(l[11])))
		
		if int(l[2]) in hit_snp_positions:
			hit_snps.append(dict(snp))
		else:
			nonhit_snps.append(dict(snp))

	file_snps.close()
	return [hit_snps, nonhit_snps]

def plotRDistribution_hit_androniki ():
	hit_r = [] 
	[hit_snps, nonhit_snps] = readSNPData_Androniki()
	for snp in hit_snps:
		hit_r.append(snp['r'])

	bins = np.linspace(-1.0,1.0,20)
	w = np.ones_like(hit_r) / len(hit_r)
	plt.hist(hit_r, bins, weights=w, facecolor='blue', alpha=1)
	plt.xlabel(r'$R$')
	plt.ylabel('Probability')
	plt.title(r'GWAS SNP - Deletion combinations on chromosome 20')
	plt.savefig('R_chr20.png')
	plt.show() 
	

def plotRDistribution_hit ():
	hit_r = [] 
	for chromosome in range(1,23):
		print 'Processing chromosome', chromosome
		[hit_snps, nonhit_snps] = readSNPData(chromosome)
		for snp in hit_snps:
			hit_r.append(snp['r'])
			#print snp['r']
	bins = np.linspace(-1.0,1.0,100)
	w = np.ones_like(hit_r) / len(hit_r)
	plt.hist(hit_r, bins, weights=w, facecolor='blue', alpha=1)
	plt.xlabel(r'$R$')
	plt.ylabel('Probability')
	plt.title(r'GWAS SNP - Deletion combinations')
	plt.savefig('R.png')
	plt.show() 


def plotGamma_n ():
	hit, nonhit = [], [] 
	for chromosome in range(22,23):
		print 'Processing chromosome', chromosome 
		[hit_snps, nonhit_snps] = readSNPData(chromosome)
		for snp in hit_snps:
			hit.append(snp['g_n'])
		for snp in nonhit_snps:
			nonhit.append(snp['g_n'])
	bins = np.linspace(-1.0,1.0,100)
	w = np.ones_like(hit) / len(hit)
	plot1 = plt.hist(hit, bins, weights=w, facecolor='red', alpha=0.5, label='GWAS SNPs')
	w = np.ones_like(nonhit) / len(nonhit)
	plot2 = plt.hist(nonhit, bins, weights=w, facecolor='blue', alpha=0.5, label='non GWAS SNPs')
	plt.xlabel(r'$\gamma_n$')
	plt.ylabel('Probability')
	plt.title('SNP-Deletion Association')
	plt.legend()
	plt.show() 


def plotMaxR2 ():
	hit, nonhit = [], [] 
	for chromosome in range(1,23):
		print 'Processing chromosome', chromosome 
		[hit_snps, nonhit_snps] = readSNPData(chromosome)
		for snp in hit_snps:
			hit.append(snp['r']**2)
		for snp in nonhit_snps:
			nonhit.append(snp['r']**2)
	bins = np.linspace(0.0,1.0,100)
	w = np.ones_like(hit) / len(hit)
	plot1 = plt.hist(hit, bins, weights=w, facecolor='red', alpha=0.5, label='GWAS SNPs')
	w = np.ones_like(nonhit) / len(nonhit)
	plot2 = plt.hist(nonhit, bins, weights=w, facecolor='blue', alpha=0.5, label='non GWAS SNPs')
	plt.xlabel(r'$R^2$')
	plt.ylabel('Probability')
	plt.title('SNP-Deletion Association')
	plt.legend()
	plt.show() 

def plotECDF (): 
	hit, nonhit = [], [] 
	for chromosome in range(1,23):
		print 'Processing chromosome', chromosome 
		snps = readSNPData(chromosome)
		for snp in snps:
			if snp['hit'] == 1:
				hit.append(snp['r2'])
			else:
				nonhit.append(snp['r2'])
	num_bins = 100
	counts, bin_edges = np.histogram(nonhit, bins=num_bins, normed=True)
	cdf = np.cumsum(counts)
	plt.plot(bin_edges[1:], cdf)
 	plt.show()



# Writes arrays to a prescripted output file 
def writeSNPsToFile(snps, filename):
	output_file = open(filename, 'a')
	output_file.write('pos\tsnp_af\ttype\tR\tp\tdel_pos\tdel_length\tdel_af\td_tss\n')
	for snp in snps:
		output_file.write(str(snp['pos']) + '\t')
		output_file.write(str(snp['af']) + '\t')	
		output_file.write(str(snp['type']) + '\t')
		output_file.write(str(snp['r']) + '\t')	
		output_file.write(str(snp['p']) + '\t')
		output_file.write(str(snp['del_pos']) + '\t')
		output_file.write(str(snp['del_length']) + '\t')
		output_file.write(str(snp['del_af']) + '\t')
		output_file.write(str(snp['d_tss']) + '\n')				

	output_file.close() 

#plotRDistribution_nonhit()
#plotMaxR2() 
#sample3()

#[hit_snps, nonhit_snps] = readSNPData(20)
#writeSNPsToFile(hit_snps, 'hitSNPsChromosome20.txt')

#for c in range(1,22):
#	compute(c, ('children'))












# Read in the data for one chromosome and returns the data for the hit- and nonhit-SNPs and deletion pairs
def readData(chromosome, individuals):
	file_snps = None
	if 'children' in individuals:
		file_snps = open('chr' + str(chromosome) + '_children_last.txt', 'r')
	else:
		file_snps = open('chr' + str(chromosome) + '_parents_last.txt', 'r')

	hit_snps, nonhit_snps = [], [] 

	# Read through all lines
	for line in file_snps.readlines():
		l = [value for value in line.split()]		
		snp = [('pos', int(l[1])), ('del_pos', int(l[2])), ('type', int(l[3])), ('d_tss', int(l[4])), ('d_del', int(l[5])), ('r', float(l[6])), ('p', float(l[7])), ('a', int(l[8])), ('b', int(l[9])), ('c', int(l[10]))]
		if int(l[0]) == 1:
			hit_snps.append(dict(snp))
		else:
			nonhit_snps.append(dict(snp))
	file_snps.close()
	return [hit_snps, nonhit_snps]

# Determines whether two SNP-Deletion pairs are similar
def snpsSimilarP (hit_snp, nonhit_snp, d_tss_thres = D_TSS_THRES, d_snp_del_thres = D_SNP_DEL_THRES):
	if hit_snp['type'] == nonhit_snp['type'] and abs(hit_snp['d_tss'] - nonhit_snp['d_tss']) <= d_tss_thres and abs(hit_snp['d_del'] - nonhit_snp['d_del']) <= d_snp_del_thres:	
		af = (hit_snp['a'] + hit_snp['b']) / 470.0 # TODO chance for parents!
		n_af = (nonhit_snp['a'] + nonhit_snp['b']) / 470.0		
		if abs(af - n_af) <= 0.05 or abs(af - (1.0 - n_af)) <= 0.05:		
			return True
		else:
			return False 
	else:
		return False 

# Returns a list that contains for every hit snp a list of indices of nonhit-snps that are similar 
def returnControlIndicesP (hit_snps, nonhit_snps, d_tss_thres = D_TSS_THRES, d_snp_del_thres = D_SNP_DEL_THRES):
	control_indices = [] # will contain lists of indices. Each list contains indices for potential control snps
	n_no_similar_nonhitsnps = 0 	

	for hit_snp in hit_snps:
		similar_nonhitsnps = [] 
		for i in range(len(nonhit_snps)):
			if snpsSimilarP(hit_snp, nonhit_snps[i], d_tss_thres, d_snp_del_thres):
				similar_nonhitsnps.append(i)
		
		if len(similar_nonhitsnps) == 0:
			n_no_similar_nonhitsnps = n_no_similar_nonhitsnps + 1
		
		control_indices.append(similar_nonhitsnps)

	print '# of hit SNPs with no similar non-hit SNPs:', n_no_similar_nonhitsnps
	return control_indices 


def fdr (p_values, q, m):
	p_values.sort()
	for k in range(1,m+1):
		if p_values[k-1] > (k / float(m)) * q:
			return k-1 
	return m 


def sampleP_final(chromosomes, d_tss_thres = D_TSS_THRES, d_snp_del_thres = D_SNP_DEL_THRES, n_samples = N_SAMPLES, individuals=('children')):
	sample_p = [[]] * n_samples # p-values are stored here for every sample
	hit_p = [] # p-values for hit snps are stored here
	n_hit_snps = 0.0 

	for chromosome in chromosomes:
		print 'Chromosome', chromosome
		[hit_snps, nonhit_snps] = readData(chromosome, individuals)
		n_hit_snps = n_hit_snps + len(hit_snps)
		control_indices = returnControlIndicesP (hit_snps, nonhit_snps, d_tss_thres, d_snp_del_thres)
		
		for hit_snp in hit_snps:
			hit_p.append(hit_snp['p'])

		for sample in range(n_samples):
			for i in range(len(hit_snps)):
				if len(control_indices[i]) != 0:
					p = nonhit_snps[choice(control_indices[i])]['p']
					sample_p[sample] = sample_p[sample] + [p]
				else:
				 	p = choice(nonhit_snps)['p']
					sample_p[sample] = sample_p[sample] + [p]

	for sample in range(n_samples):
		sample_p[sample] = fdr(sample_p[sample], .05, int(n_hit_snps))

	case = fdr(hit_p, .05, int(n_hit_snps))
	return [n_hit_snps, case, sample_p]

def sampleP (chromosome, d_tss_thres = D_TSS_THRES, d_snp_del_thres = D_SNP_DEL_THRES, n_samples = N_SAMPLES, individuals=('children')):
	[hit_snps, nonhit_snps] = readData(chromosome, individuals)
	n_hit_snps = len(hit_snps)
	print n_hit_snps
	# control_indices contain lists of indices. Each list contains indices for potential control snps
	control_indices = returnControlIndicesP (hit_snps, nonhit_snps, d_tss_thres, d_snp_del_thres)
	print len(control_indices)
	sample_results = [None] * n_samples 

	for sample in range(n_samples):
		p_values = [] 
		for i in range(n_hit_snps):
			if len(control_indices[i]) != 0:
				p_values.append(nonhit_snps[choice(control_indices[i])]['p'])
			else:
				p_values.append(choice(nonhit_snps)['p'])
		sample_results[sample] = fdr(p_values, .10, n_hit_snps)		
		print sample_results[sample] 	

	p_values_hit = []
	for hit_snp in hit_snps:
		p_values_hit.append(hit_snp['p'])
	case = fdr(p_values_hit, .10, n_hit_snps)	
	return [n_hit_snps, case, sample_results]

print sampleP_final(range(1,23), d_tss_thres = 10000, d_snp_del_thres = 25000, n_samples = 10000, individuals=('parents'))

