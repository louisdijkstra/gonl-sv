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

import vcf as vcf # working with vcf files 
import math 
from scipy.stats.stats import pearsonr
import scipy as sp
import scipy.stats


DIST_THRESHOLD = 1000000 # 1M bp distance threshold between deletions and SNPs
MAJOR_AF_THRESHOLD = .96 # no restriction on which deletions to consider
P_THRESHOLD = 0.05 
R2_THRESHOLD = 0.8

# Returns the distance between the SNP and the deletion (CNV)
# Returns 0 when the SNP lies within the deletion
def distanceSNP_DEL (snp_pos, del_start, del_length):
	if snp_pos < del_start:
		return del_start - snp_pos
	elif snp_pos > del_start + del_length:
		return snp_pos - del_start + del_length 
	else: # snp overlaps with the deletion
		return 0 

# determine Pearson's R for a 2x2 table of the form: 
#	a	  c 	| a + c
#	b 	  d 	| b + d
#	-------------------------
#	a + b 	c + d	| a+b+c+d
def determineR(a,b,c,d):
	return float((a * d - b * c)) / math.sqrt((a+b)*(a+c)*(b+d)*(c+d))

# Returns phase data and the major allele frequency 
def determinePhase(line, indices, hitallele):
	phase = [] 	
	for i in indices:
		phase.append(int(line[i].split('|')[0]))
		phase.append(int(line[i].split('|')[1]))
	
	if line[3] == hitallele:
		for i in range(len(phase)):
			phase[i] = 1 - phase[i]
	
	maf = sum(phase) / float(len(phase))
	if maf < .5:
		return [phase, 1 - maf]
	else:
		return [phase, maf]	


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
	

# obtain SNP data for all snps with the positions given in 'snp_list'
def obtainSNPData (vcf_filename, snp_list, indices_individuals):
	vcf_file = vcf.openVCFFile(vcf_filename)
	vcf.discardVCFHeaders(vcf_file)
	vcf_file.readline() # discard the column description header

	snps = [] # data is stored here
	for variant in vcf_file.readlines(): 
		data = [value for value in variant.split()]
		snp_pos = int(data[1])
		if len(data[3]) == 1 and snp_pos in snp_list: # variant is a SNP & present in the array
			[phase, maf] = determinePhase (data, indices_individuals, None)
			if maf <= MAJOR_AF_THRESHOLD: 
				snps.append(dict([('pos', snp_pos), ('phase', phase), ('maf', maf)]))

	vcf_file.close()
	return snps

# returns the data necessary for determining the tag SNPs for each deletion
def obtainData (chromosome, individuals):
	# create the vcf file object
	vcf_filename = returnVCFFileName (chromosome)
	# retrieve indices of the individuals 
	indices_individuals = vcf.returnColumns (vcf_filename, individuals) 
	n_individuals = len(indices_individuals)
	# obtain the deletions found on the chromosome 
	deletions = vcf.returnDeletions (chromosome, indices_individuals, MAJOR_AF_THRESHOLD)
	# obtain the list of SNPs that are to be considered (the ones on a specific array)
	snp_list = obtainSNPList(chromosome)
	# obtain relevant SNP data from vcf file (only for the SNPs in snp_list)
	snps = obtainSNPData (vcf_filename, snp_list, indices_individuals)
	return [deletions, snps, indices_individuals, n_individuals]


# Writes arrays to a prescripted output file 
def writeToOutputFile(output, individuals):
	filename = None
	if 'children' in individuals:
		filename = 'tagSNPs_for_deletions_children.txt' 
	else:
		filename = 'tagSNPs_for_deletions_parents.txt'  
	output_file = open(filename, 'a')
	for line in output:
		for item in line:
			output_file.write(str(item) + '\t')
		output_file.write('\n')

	output_file.close() 


def returnVCFFileName (chromosome):
	if chromosome == 1 or chromosome == 22: 
		return '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chromosome) + '_snps_and_dels.vcf'
	else: 
		return '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chromosome) + '_snps_and_dels.vcf.gz'

def determineTagSNPsChromosome (chromosome, individuals=('parents')):
	# collect the necessary data
	[deletions, snps, indices_individuals, n_individuals] = obtainData (chromosome, individuals)
	
	tagSNPs = [] # tag snp data is stored here
	
	# walk through the deletions 
	for deletion in deletions: 
		candidate_snps = [] # list of candidate tag-snps
		# walk through all snps
		for snp in snps:
			# determine whether the distance between SNP and DEL is below the preset threshold
			dist_snp_del = distanceSNP_DEL (snp['pos'], deletion['pos'], deletion['length'])
			if dist_snp_del <= DIST_THRESHOLD:
				# determine the 2x2 table and compute Fisher's p and R
				[a,b,c,d] = vcf.return2x2Table(snp['phase'], deletion['phase'])
				[odds_ratio, p] = sp.stats.fisher_exact([[a,c],[b,d]], alternative='two-sided')
				r = determineR(a,b,c,d)
				prob_a = a / float(a + b) 
				prob_not_a = c / float(c + d) 
				candidate_snps.append([snp['pos'], r, p, prob_a, prob_not_a])
		
		tagSNP = None 
		bestR2 = R2_THRESHOLD
		# walk through the candidate snps and select the 'best' one
		for candidate in candidate_snps: 
			if candidate[2] <= P_THRESHOLD and candidate[1]**2 >= bestR2: 
				tagSNP = candidate 
		
		if tagSNP is not None:
			tagSNPs.append([chromosome, deletion['pos'], deletion['length']] + tagSNP)
			print [chromosome, deletion['pos'], deletion['length']] + tagSNP
		else: 
			print 'No suitable tag SNPs were found for deletion on position ', deletion['pos']
		
	writeToOutputFile(tagSNPs, individuals)

# determine the tag SNPs for all autosomes
def determineTagSNPs (individuals=('parents')):
	for chromosome in range(1,22):
		determineTagSNPsChromosome (chromosome, individuals)


#determineTagSNPsChromosome (22, individuals=('parents'))
determineTagSNPs(individuals=('parents'))



		
