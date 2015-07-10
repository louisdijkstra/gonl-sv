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

# processPhasedData.py processes the phased data of SNPs and deletions
#	from the GoNL project. Here, we focus only on the data of 
#	independent children (we consider only one child from each twin pair).
#	
#	For every pair (deletion, SNP) which are less than 2 million bp 
#	away from each other, we compute Pearson's correlation (r) and 
#	the corresponding p-value. 
#
#	For every autosome (1-22) the output is stored in a file of the 
#	form: 
#
#		del_snp_pearson_r_chr?.phased
#
#	Each file consists of sections, where each section contains the 
# 	data associated with one deletion. The section starts with a line 
#	enumerating the 7 data points associated with the deletion: 
#
#		#CHR POS LENGTH MAF 	
#
#	where 
#		#CHR - the chromosome where the deletion took place
#		POS - the start position of the deletion
#		LENGTH - the estimated length of the observed deletion
#		MAF - major allele frequency
#
#	A list of SNPs that lie within a radius of 2M bp of the deletion
#	follows the 'deletion-line':
#
#		#CHR POS TYPE MAF DIST_DEL DIST_TSS R P 
#
#	where 
#		#CHR - the chromosome on which the SNP lies
#		POS - the position of the SNP
#		TYPE - intergenic/intronic/exonic
#		MAF - major allele frequency
#		DIST_DEL - distance to the deletion
#		DIST_TSS - minimal distance to TSS
#		R - Pearson's correlation coefficient
#		P - corresponding p-value 
# 
#	Every deletion section is ended with a line with only a dot '.'
# 
# Author: Louis Dijkstra
# Date: October 2013
import math 
import gzip
from scipy.stats.stats import pearsonr
import sys
from time import sleep

DIST_THRES = 2000000 


# Returns the distance between the SNP and the deletion (CNV)
# Returns -1 and prints a warning when the SNP lies within the deletion
def distanceSNP_CNV (snp_pos, cnv_start, cnv_length):
	if snp_pos < cnv_start:
		return cnv_start - snp_pos - 1
	elif snp_pos > cnv_start + cnv_length - 1:
		return snp_pos - cnv_start - cnv_length 
	else:
		print 'WARNING: SNP at position ', snp_pos, ' appears to have been deleted' 
		return -1 

# Returns the type of the SNP (intergenic/intronic/exonic)given 
# the chromosomal data (chromosome)
def determineType (snp_pos, chromosome):
	# move through all coding regions in the chromosome
	for gene in chromosome:
		if gene[2] <= snp_pos <= gene[3]: # SNP is in coding region
			# determine whether SNP is in exon 
			for exon in range(gene[4]):
				# whether SNP lies within the exon
				if gene[5][exon][0] <= snp_pos <= gene[5][exon][1]: 
					return 'exonic'
			return 'intronic'

	return 'intergenic'

# Returns the distance to the nearest TSS given chromosomal data (chromosome)
def distanceToTSS (snp_pos, chromosome):
	min_dist = float('inf') 

	for gene in chromosome:
		if abs(gene[0] - snp_pos) < min_dist:
			min_dist = abs(gene[0] - snp_pos)  

	return min_dist 

# Returns the major allele frequency (MAF) 
def determineMAF (af):
	if af >= .5: 
		return af
	else:
		return 1.0 - af


# Prints chromosomal data to the screen 
def printChromosome (chromosome):
	print 'TSS_start\tTSS_end\tcdsStart\tcdsEnd\tn_exon\texonStarts\texonEnds\n'
	for gene in chromosome:
		for i in gene:
			print i, '\t', 
		print '\n'

# Writes one array to predescripted output file
def writeLineToOutputFile(output, chr):
	filename = 'del_snp_pearson_r_chr' + str(chr) + '.phased' 
	output_file = open(filename, 'a')
	for item in output:
		output_file.write(str(item) + '\t')
	output_file.write('\n')

	output_file.close() 

# Writes arrays to a prescripted output file 
def writeToOutputFile(output, chr):
	filename = 'del_snp_pearson_r_chr' + str(chr) + '.phased' 
	output_file = open(filename, 'a')
	for line in output:
		for item in line:
			output_file.write(str(item) + '\t')
		output_file.write('\n')

	output_file.close() 

# Reads chromosomal data from file for a specified chromosome (chr)
# Each line represents one gene. Columns represent:
#	1) transcription start site (TSS)
#	2) transcription end site 
#	3) start coding region
#	4) end coding region
#	5) number of exons
#	6) list with exons starts and exon ends
def readChromosomeData (chr):

	chr = str(chr)		
	print 'Reading data for chromosome', chr, '...'
	
	hg19_file = open('/ufs/dijkstra/Projects/SNPs_LD_deletions_phased/data/hg19.txt', 'r') # open data file
	hg19_file.readline() # discard the header 
	
	chromosome = [] ; # data is stored here

	# Read through all lines...
	for line in hg19_file.readlines():
		array = [value for value in line.split()]
		if array[2] == 'chr' + chr:
			gene = map(int, array[4:9])
			exon_start = array[9].split(',')
			exon_end = array[10].split(',')
			gene.append(zip(map(int, exon_start[:gene[-1]]), map(int, exon_end[:gene[-1]])))
			chromosome.append(gene)

	hg19_file.close()			
	return chromosome

# Returns the (zero-counting) index of the columns related to children 
def retrieveColumnsOfChildren (chr):
	# open the file with phased data
	phased_file = gzip.open('/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chr) + '_snps_and_dels.vcf.gz', 'r')

	# discard the headers (first 94 lines)
	for line in range(94):
		phased_file.readline() 

	columns_children = [] 

	header = [value for value in phased_file.readline().split()] 

	for i in range(len(header)):
		if header[i][-1] == 'c':
			columns_children.append(i)

	phased_file.close() 
	print 'Number of children chr' + str(chr) + ':',  len(columns_children)
	return columns_children


# Retrieves the deletions for a file:
# Returns:
#	deletions - every line a deletion
#		POS
#		LENGTH
#		Major allele frequency
#		Vector
def retrieveDeletions (chr, child):

	# open the file with phased data
	phased_file = gzip.open('/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chr) + '_snps_and_dels.vcf.gz', 'r')

	# discard the headers (first 94 lines)
	for line in range(95):
		phased_file.readline() 

	n_dels = 0 
	n_dels_with_maf1 = 0 

	deletions = [] # data is stored here

	# loop over all lines until you find a deletion
	for line in phased_file.readlines():
		l = [value for value in line.split()]
		if len(l[3]) > 1: # deletion! 
			# pos / length 
			deletion = [int(l[1]), len(l[3]) - 1]
			# retrieve AF
			del_maf = determineMAF(float(l[7].split(';')[1].split('=')[1]))
			deletion.append(del_maf)
			
			# create vector:
			x = [] 
			for c in child:
				k = l[c].split('|')
				x.append(float(k[0]))
				x.append(float(k[1]))
			
			deletion.append(x)
			if del_maf < 1: 
				deletions.append(deletion)
				n_dels = n_dels + 1
			else:
				n_dels_with_maf1 = n_dels_with_maf1 + 1


	phased_file.close()
	print '---------------------------------' 
	print 'Chromosome ' + str(chr) + ': '
	print 'Deletions with major allele frequency < 1: ', n_dels
	print 'Deletions with major allele frequency = 1: ', n_dels_with_maf1
	print 'Total number of deletions: ', n_dels + n_dels_with_maf1
	print '---------------------------------' 
	return deletions 

def processPhasedData (chr, child, chromosome):

	dels = retrieveDeletions (chr, child) 

	# open the file with phased data
	phased_file = gzip.open('/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chr) + '_snps_and_dels.vcf.gz', 'r')
	output = [] # results for this deletion are stored here
	# discard the headers (first 94 lines)
	for line in range(95):
		phased_file.readline() 

	iSNP = 0 
	for line in phased_file.readlines():
		l = [value for value in line.split()]
		# check whether line contains SNP data
		snp_maf = determineMAF(float(l[7].split(';')[1].split('=')[1]))
		if snp_maf <= .96 and len(l[3]) == 1:
			snp_pos = int(l[1]) 
			snp_type = determineType(snp_pos, chromosome)				
			dist_tss = distanceToTSS (snp_pos, chromosome)
			# create vector:
			y = [] 
			for c in child:
				k = l[c].split('|')
				y.append(float(k[0]))
				y.append(float(k[1]))
			
			pos_del = 0 
			length_del = 0 
			maf_del = 0 
			dist_del = 0 
			maxR2 = 0 
			P_final = 0 

			for d in dels: 
				if abs(snp_pos - d[0]) <= DIST_THRES: 
					[R2, P] = pearsonr(d[3], y) 	
					R2 = R2 ** 2 		
					if not math.isnan(R2) and maxR2 < R2: 
						maxR2 = R2 
						P_final = P 
						pos_del = d[0]
						length_del = d[1]
						maf_del = d[2]
						dist_del = distanceSNP_CNV(snp_pos, pos_del, length_del)
			
			output.append([chr, snp_pos, snp_type, snp_maf, pos_del, length_del, maf_del, maxR2, P_final, dist_del, dist_tss]) 
			iSNP = iSNP + 1 			
	
	print 'Chromosome ' + str(chr) + '| number of SNPs: ', iSNP 
	writeToOutputFile(output, chr) 
	phased_file.close() 

def processPhasedData2 (chr, child, chromosome):
	
	dels = retrieveDeletions (chr, child) 
	for d in dels:
		del_pos = d[0]
		del_len = d[1]
		del_vec = d[3] 

		print [chr, del_pos, del_len, d[2]]

		writeLineToOutputFile([chr, del_pos, del_len, d[2]], chr)

		# open the file with phased data
		phased_file = gzip.open('/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chr) + '_snps_and_dels.vcf.gz', 'r')
		output = [] # results for this deletion are stored here
		# discard the headers (first 94 lines)
		for line in range(95):
			phased_file.readline() 
			
		for line in phased_file.readlines():
			l = [value for value in line.split()]
			# SNP within the preset range 
			if len(l[3]) == 1 and abs(int(l[1]) - del_pos) <= DIST_THRES:
				snp_pos = int(l[1]) 
				snp_type = determineType(snp_pos, chromosome)				
				snp_maf = determineMAF(float(l[7].split(';')[1].split('=')[1]))
				dist_del = distanceSNP_CNV (snp_pos, del_pos, del_len) 
				dist_tss = distanceToTSS (snp_pos, chromosome)
				# create vector:
				y = [] 
				for c in child:
					k = l[c].split('|')
					y.append(float(k[0]))
					y.append(float(k[1]))
				[R,P] = pearsonr(del_vec, y) 
				if not math.isnan(R):
					output.append([chr, snp_pos, snp_type, snp_maf, dist_del, dist_tss, R, P])  
				
		print output[-1]
		writeToOutputFile(output, chr) 
		phased_file.close() 
		writeLineToOutputFile(['---'], chr)
	
def main():
	for chr in range(2,22):
		child = retrieveColumnsOfChildren(chr)
		chromosome = readChromosomeData(chr)
		processPhasedData(chr, child, chromosome)

main()



