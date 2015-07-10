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

# File to obtain the necessary data for the paper
# Full of small functions that incorporates data from several files in /data

import math 
from scipy.stats.stats import pearsonr
import vcf as vcf 
import hg19 as hg19
import gamma as gamma


# reads the hit SNPs from GWAS_catalogue_caucasian_SNPIDs.alleles.txt
# Hit SNPs for which the hit allele is not defined are discarded and stored in discarded_hit_snps_positions
def readHitSNPs (chromosome):
	snp_files = open('/ufs/dijkstra/Projects/SNPs_CNVs/data/GWAS_catalogue_caucasian_SNPIDs.alleles.txt', 'r')
	snp_files.readline() # discard header
	
	hit_snps_positions = [] # positions of the hit snps are stored here
	snps = [] # positions and hit allele are stored here
	discarded_hit_snps_positions = [] 

	chromosome = 'chr' + str(chromosome)

	for line in snp_files.readlines():
		info = [value for value in line.split()]
		if info[0] == chromosome: # snp lies on the chromosome of interest
			if len(info) == 4: # check whether the hit allele is given in the file
				hit_snps_positions.append(int(info[1]))
				snps.append( dict( [('pos', int(info[1])), ('hit_allele', info[3])] ) )
			else:
				discarded_hit_snps_positions.append(int(info[1]))
	
	snp_files.close()
	return [hit_snps_positions, snps, discarded_hit_snps_positions]

# returns the list of snps that are to be considered, i.e., that are on the array of interest
def obtainSNPList (chromosome):
	chromosome = str(chromosome)
	HumanOmni_file = open('/ufs/dijkstra/Projects/SNPs_CNVs/data/HumanOmni25-8v1-1_A-b37.strand', 'r')
	snp_list = [] 
	for line in HumanOmni_file.readlines():
		l = [value for value in line.split()]
		if l[1] == chromosome:
			snp_list.append(int(l[2]))
	return snp_list



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

# return number of SNPs on the HumanOmni25 list. Only autosomal SNPs!
def returnNumberOfSNPsOnArray ():
	number_of_snps = 0 
	for chromosome in range(1,23):
		snp_list = obtainSNPList (chromosome)
		number_of_snps = number_of_snps + len(snp_list)

	print 'Number of autosomal SNPs on HumanOmni25', number_of_snps
	return number_of_snps

# return number of hit SNPs considered 
def returnNumberOfHitSNPs ():
	n = 0 
	for chromosome in range(1,23):
		[hit_snps_positions, snps, discarded_hit_snps_positions] = readHitSNPs (chromosome)
		n = n + len(hit_snps_positions)
	print 'Number of hit SNPs originally considered', n
	return n


def returnNumberOfSNPsDeletionsMAFChr (chromosome):
	file_snps = open('chr' + str(chromosome) + '_parents_last.txt', 'r')
	
	list_hsnps = [] 
	list_nhsnps = [] 
	list_deletions = [] 
	for line in file_snps.readlines():
		l = [value for value in line.split()]	
		if int(l[0]) == 1: # hit snp
			list_hsnps.append(int(l[1]))
		else:
			list_nhsnps.append(int(l[1]))
		list_deletions.append(int(l[2]))
	file_snps.close()
	list_hsnps = set(list_hsnps)
	list_nhsnps = set(list_nhsnps)
	list_deletions = set(list_deletions)
	print chromosome, len(list_hsnps), len(list_nhsnps),len(list_deletions)
	return [len(list_hsnps), len(list_nhsnps), len(list_deletions)]
	
def returnNumberOfSNPsDeletionsWithMAF ():
	n_hsnps = 0 
	n_nhsnps = 0 
	n_deletions = 0 
	for chromosome in range(1,23):
		[n_s, n_nhs, n_d] = returnNumberOfSNPsDeletionsMAFChr (chromosome)
		n_hsnps = n_hsnps + n_s
		n_nhsnps = n_nhsnps + n_nhs
		n_deletions = n_deletions + n_d

	print 'Number of SNPs and deletions after filtering out everything with MAF < .04', n_hsnps, n_nhsnps, n_deletions


returnNumberOfSNPsDeletionsWithMAF ()
	




