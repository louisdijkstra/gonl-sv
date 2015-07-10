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
import vcf as vcf 
import hg19 as hg19

from operator import itemgetter


# Global variables
LD_THRES = .8 	# R^2 threshold for a deletion and SNP to be considered 'related'
D_TSS_THRES = 1000 
D_SNP_DEL_THRES = 100000
AF_THRES = 0.05
N_SAMPLES = 1000 # number of samples 




# Read in the data for one chromosome and returns the data for the hit- and nonhit-SNPs and deletion pairs
def readData(chromosome, individuals):
	file_snps = None
	if 'children' in individuals:
		file_snps = open('chr' + str(chromosome) + '_children_last.txt', 'r')
	else:
		file_snps = open('chr' + str(chromosome) + '_parents_last.txt', 'r')

	hit_snps = [] 

	# Read through all lines
	for line in file_snps.readlines():
		l = [value for value in line.split()]	
		if int(l[0]) == 1: # hit SNP
			snp = [chromosome, int(l[1]), ]
			snp = [('pos', int(l[1])), ('del_pos', int(l[2])), ('type', int(l[3])), ('d_tss', int(l[4])), ('d_del', int(l[5])), ('r', float(l[6])), ('p', float(l[7])), ('a', int(l[8])), ('b', int(l[9])), ('c', int(l[10]))]
		if int(l[0]) == 1:
			hit_snps.append(dict(snp))
		else:
			nonhit_snps.append(dict(snp))
	file_snps.close()
	return [hit_snps, nonhit_snps]




	

# reads the hit SNPs from GWAS_catalogue_caucasian_SNPIDs.alleles.txt
# Hit SNPs for which the hit allele is not defined are discarded and stored in discarded_hit_snps_positions
def readHitSNPData (chromosome):
	snp_files = open('/ufs/dijkstra/Projects/SNPs_CNVs/data/GWAS_catalogue_caucasian_SNPIDs.alleles.txt', 'r')
	snp_files.readline() # discard header
	
	hit_snps_positions = [] # positions of the hit snps are stored here
	snps = [] # positions and hit allele are stored here

	chromosome = 'chr' + str(chromosome)

	for line in snp_files.readlines():
		info = [value for value in line.split()]
		if info[0] == chromosome: # snp lies on the chromosome of interest
			if len(info) == 4: # check whether the hit allele is given in the file
				hit_snps_positions.append(int(info[1]))
				snps.append( dict( [('pos', int(info[1])), ('hit_allele', info[3])] ) )

	
	snp_files.close()
	return [hit_snps_positions, snps]


def returnVCFFileName (chromosome):
	if chromosome == 1 or chromosome == 22: 
		return '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chromosome) + '_snps_and_dels.vcf'
	else:
		return '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chromosome) + '_snps_and_dels.vcf.gz'



def retrieveAllele (snp_pos, allele_info):
	hit_allele = None
	for item in allele_info:
		if item['pos'] == snp_pos:
			hit_allele = item['hit_allele'] 
	if hit_allele == None:
		print 'ERROR: hit allele could not be found for some reason....'
	return hit_allele 

def retrieveDelLength(del_pos, deletions):
	length = None 
	for item in deletions:
		if item['pos'] == del_pos:
			length = item['length']
	if length == None:
		print 'ERROR: deletion could not be found for some reason...'
	return length

def collectAllData (chromosome, individuals, deletions, allele_info, n_individuals):
	file_snps = None
	if 'children' in individuals:
		file_snps = open('chr' + str(chromosome) + '_children_last.txt', 'r')
	else:
		file_snps = open('chr' + str(chromosome) + '_parents_last.txt', 'r')

	hit_snps = [] 

	for line in file_snps.readlines():
		l = [value for value in line.split()]	
		if int(l[0]) == 1: # hit SNP
			snp_pos = int(l[1])
			hit_allele = retrieveAllele(snp_pos, allele_info)
			del_pos = int(l[2])
			del_length = retrieveDelLength(del_pos, deletions)
			dist_tss = int(l[4])
			r = float(l[6])			
			p = float(l[7])
			a, b, c = int(l[8]), int(l[9]), int(l[10])
			d = 2 * n_individuals - a - b - c
			hit_snps.append([chromosome, snp_pos, hit_allele, del_pos, del_length, dist_tss, r, p, a, b, c, d])
			
	return hit_snps



def fdr(hit_snps, q):
	hit_snps.sort(key=itemgetter(7)) # sort the hit snps in ascending order for the p-values
	m = len(hit_snps) # number of hypotheses
	sign_hit_snps = [] 
	for k in range(m):
		p = hit_snps[k][7]
		if p <= (float(k + 1) / float(m)) * q:
			sign_hit_snps.append(hit_snps[k])
	return sign_hit_snps

def materialSignificanceTest (hit_snps):
	sign_hit_snps = []
	for snp in hit_snps:
		if snp[6]**2 >= .8:
			sign_hit_snps.append(snp)
	return sign_hit_snps 

# Writes arrays to a prescripted output file 
def writeToOutputFile(output, individuals):
	filename = None
	if 'children' in individuals:
		filename = 'significant_hitSNPs_children.txt' 
	else:
		filename = 'significant_hitSNPs_parents.txt'  
	output_file = open(filename, 'a')
	for line in output:
		for item in line:
			output_file.write(str(item) + '\t')
		output_file.write('\n')

	output_file.close() 

def obtainSignificantHitSNPs (chromosomes, individuals=('parents')):
	hit_snps = [] # list of all significant hit snps 

	# process each chromosome in the list chromosomes
	for chromosome in chromosomes:
		print 'Processing chromosome ', chromosome
		# obtain the SNP data (position and hit allele)
		[hit_snps_positions, allele_info] = readHitSNPData (chromosome)
		# create the vcf file object
		vcf_filename = returnVCFFileName (chromosome)
		# retrieve indices of the individuals 
		indices_individuals = vcf.returnColumns (vcf_filename, individuals) 
		n_individuals = len(indices_individuals)		
		# retrieve the deletions on the chromosome 
		deletions = vcf.returnDeletions (chromosome, indices_individuals, .96) 
		# retrieve the already processed hit snp data and incorporate the relevant info 
		hit_snps_chr = collectAllData (chromosome, individuals, deletions, allele_info, n_individuals)
		# add the snps for this chromosome to the list 
		for hit_snp in hit_snps_chr:
			hit_snps.append(hit_snp)

	print 'Before FDR: ', len(hit_snps)
	# apply fdr 
	hit_snps = fdr(hit_snps, 0.05)
	print 'After FDR: ', len(hit_snps)
		
	# apply material significance level (r^2 >= .8)
	hit_snps = materialSignificanceTest (hit_snps)
	print 'After Testing for material significance (R^2 >= .8): ', len(hit_snps)
	writeToOutputFile (hit_snps, individuals)

obtainSignificantHitSNPs (range(1,23), individuals=('children'))

