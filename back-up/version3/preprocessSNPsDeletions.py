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
from scipy.stats.stats import pearsonr
import vcf as vcf 
import hg19 as hg19
import gamma as gamma
import scipy as sp
import scipy.stats

DIST_THRESHOLD 		= 1000000
MAJOR_AF_THRESHOLD 	= .96 

# Returns the distance between the SNP and the deletion (CNV)
# Returns 0 when the SNP lies within the deletion
def distanceSNP_DEL (snp_pos, del_start, del_length):
	if snp_pos < del_start:
		return del_start - snp_pos
	elif snp_pos > del_start + del_length:
		return snp_pos - del_start + del_length 
	else: # snp overlaps with the deletion
		return 0 


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

# returns the list of snps that are to be considered
def obtainSNPList (chromosome):
	chromosome = str(chromosome)
	HumanOmni_file = open('/ufs/dijkstra/Projects/SNPs_CNVs/data/HumanOmni25-8v1-1_A-b37.strand', 'r')
	snp_list = [] 
	for line in HumanOmni_file.readlines():
		l = [value for value in line.split()]
		if l[1] == chromosome:
			snp_list.append(int(l[2]))
	return snp_list
	

# reads in the necessary data for preprocessing the data
def gatherData (vcf_filename, chromosome, individuals):
	[hit_snps_positions, hit_snps, discarded_snps_positions] = readHitSNPs (chromosome) 
	genes = hg19.read (chromosome)
	indices_individuals = vcf.returnColumns (vcf_filename, individuals) 
	deletions = vcf.returnDeletions (chromosome, indices_individuals, MAJOR_AF_THRESHOLD)
	snp_list = obtainSNPList(chromosome)
	return [hit_snps_positions, hit_snps, discarded_snps_positions, snp_list, genes, indices_individuals, deletions]


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


# Writes arrays to a prescripted output file 
def writeToOutputFile(output, chr, individuals):
	filename = None
	if 'children' in individuals:
		filename = 'chr' + str(chr) + '_children_last.txt' 
	else:
		filename = 'chr' + str(chr) + '_parents_last.txt'  
	output_file = open(filename, 'a')
	for line in output:
		for item in line:
			output_file.write(str(item) + '\t')
		output_file.write('\n')

	output_file.close() 


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
	

def determineR(a,b,c,d):
	return float((a * d - b * c)) / math.sqrt((a+b)*(a+c)*(b+d)*(c+d))


def preprocessChromosome(vcf_filename, chromosome, individuals=('children','parents')):
	[hit_snps_positions, hit_snps, discarded_snps_positions, snp_list, genes, indices_individuals, deletions] = gatherData(vcf_filename, chromosome, individuals)	
	
	vcf_file = vcf.openVCFFile(vcf_filename)
	vcf.discardVCFHeaders(vcf_file)
	vcf_file.readline() # discard the column description header

	output = [] 

	for variant in vcf_file.readlines():
		data = [value for value in variant.split()]
		snp_pos = int(data[1])
		if len(data[3]) == 1 and snp_pos in snp_list and snp_pos not in discarded_snps_positions: # variant is a SNP and no deletion
			phase = [] 		
			maf = None 	
			hit = 0 
			if snp_pos in hit_snps_positions:
				hit = 1
				for i in range(len(hit_snps)):
					if hit_snps[i]['pos'] == snp_pos:			
						[phase, maf] = determinePhase(data, indices_individuals, hit_snps[i]['hit_allele'])
						#del hit_snps[i]
						break 
			else:
				[phase, maf] = determinePhase(data, indices_individuals, None)
			if maf <= MAJOR_AF_THRESHOLD:
				snp_type = hg19.determineType(snp_pos, genes)
				dist_tss = hg19.distanceToTSS(snp_pos, genes)
				for deletion in deletions:
					#if snp_pos <= deletion['pos'] + deletion['length'] + DIST_THRESHOLD: 
					dist_snp_del = distanceSNP_DEL(snp_pos, deletion['pos'], deletion['length'])
					if dist_snp_del <= DIST_THRESHOLD:
						[a,b,c,d] = vcf.return2x2Table(phase, deletion['phase'])
						[odds_ratio, p] = sp.stats.fisher_exact([[a,c],[b,d]], alternative='two-sided')
						r = determineR(a,b,c,d)
						output.append([hit, snp_pos, deletion['pos'], snp_type, dist_tss, dist_snp_del, r, p, a, b, c])
						#print [hit, snp_pos, deletion['pos'], snp_type, dist_tss, dist_snp_del, r, p, a, b, c]
					#else:
					#	break
	writeToOutputFile(output, chromosome, individuals)
	vcf_file.close()

def main():
	print 'Chromosome 1'
	vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr1_snps_and_dels.vcf'
	preprocessChromosome(vcf_filename, 1, individuals=('children'))
	for chr in range(2,22):
		print 'Chromosome ' + str(chr)
		vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chr) + '_snps_and_dels.vcf.gz'
		preprocessChromosome(vcf_filename, chr, individuals=('children'))
	#print 'Chromosome 22'
	#vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr22_snps_and_dels.vcf'
	#preprocessChromosome(vcf_filename, 22, individuals=('children'))

	print 'Chromosome 1'
	vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr1_snps_and_dels.vcf'
	preprocessChromosome(vcf_filename, 1, individuals=('parents'))
	for chr in range(2,22):
		print 'Chromosome ' + str(chr)
		vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chr) + '_snps_and_dels.vcf.gz'
		preprocessChromosome(vcf_filename, chr, individuals=('parents'))
	print 'Chromosome 22'
	vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr22_snps_and_dels.vcf'
	preprocessChromosome(vcf_filename, 22, individuals=('parents'))

#main()
vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr22_snps_and_dels.vcf'
indices_individuals = vcf.returnColumns (vcf_filename, ('children')) 
print len(indices_individuals)

indices_individuals = vcf.returnColumns (vcf_filename, ('parents')) 
print len(indices_individuals)

