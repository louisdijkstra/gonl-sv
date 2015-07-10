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

# preprocessAssociationTests.py preprocesses the results of the 
#	association tests between each SNP and every deletion 
#	less than 2M bp away from that SNP. 
#
#	For every SNP, we determine the maximum Pearson's R^2,
#	whether it is 1) introgenic, 2) intronic, or, 3) exonic,
#	its distance to the nearest TSS, its minor allele 
#	frequency (MAF), and the distance to the deletion which
#	resulted in the highest R^2. 
#
#	For every chromosome, we output a file in the form:
#
#		 SNP_max_r2_chr?.txt
#	
#	Every line represents one SNP. The columns represent:
#	
#		1) chromosome
#		2) 1 - hit/GWAS SNP, 0 - otherwise
#		3) SNP position
#		4) type (intergenic, intronic, exonic)
#		5) minor allele frequency
#		6) CNV position
#		7) CNV length
#		8) max. R^2
#		9) corresponding p-value
#		10) distance to deletion
#		11) distance to nearest TSS
# 
# Author: Louis Dijkstra (dijkstra@cwi.nl)
# Date: October 2013
import math

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

# Returns the minor allele frequency (MAF) of a SNP 
def determineSNPMAF (snp_af):
	if snp_af <= .5: 
		return snp_af
	else:
		return 1 - snp_af

# Prints chromosomal data to the screen 
def printChromosome (chromosome):
	print 'TSS_start\tTSS_end\tcdsStart\tcdsEnd\tn_exon\texonStarts\texonEnds\n'
	for gene in chromosome:
		for i in gene:
			print i, '\t', 
		print '\n'

# Prints the results of the association tests
def printResultsAssociationTests (at):
	for test in at:
		for i in test:
			print str(i) + '\t',
		print '\n'

# Writes data to a prescripted output file 
def writeToOutputFile(output, chr):
	filename = 'SNP_max_r2_chr' + str(chr) + '.txt'
	output_file = open(filename, 'a')
	for line in output:
		for item in line:
			output_file.write(str(item) + '\t')
		output_file.write('\n')

	output_file.close() 

# Reads chromosomal data from file for a specified chromosome (chr)
# Each line represents one gene. Columns represent:
#	1) transcription start site
#	2) transcription end site
#	3) start coding region
#	4) end coding region
#	5) number of exons
#	6) list with exons starts and exon ends
def readChromosomeData (chr):

	chr = str(chr)		
	print 'Reading data for chromosome', chr, '...'
	
	hg19_file = open('/ufs/dijkstra/Projects/SNPs_LD_deletions/data/hg19.txt', 'r') # open data file
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

# Reads the association tests and processes the data 
# First, we work through the hit SNPs and then process the 
# rest
def processAT(chr, chromosome):

	# read the association test results for the hit-SNPs
	at_file = open('/data2/as/hitSNP-CNV/results/HitSNP.Del-20130807.LD-2000000.causal.pearson.chr' + str(chr) + '.txt')
	at_file.readline() # discard the header
	
	output = [] # results are stored here
	
	line 		= at_file.readline().split()
	snp_pos 	= int(line[1])
	cnv_pos 	= int(line[3])
	cnv_length 	= int(line[4])
	snp_maf 	= determineSNPMAF (float(line[5]))
	max_r2 		= float(line[7])
	p 		= float(line[8])

	hit_snp_pos = [snp_pos] 

	# loop over all association test results
	for test in at_file.readlines():
		line = [value for value in test.split()]
		if snp_pos == int(line[1]): # still same SNP?
			r2 = float(line[7])
			if not math.isnan(r2) and max_r2 < r2:
				max_r2 		= r2
				cnv_pos 	= int(line[3])
				cnv_length 	= int(line[4])
				p 		= float(line[8])
		else:
			output.append([chr, 1, snp_pos, determineType(snp_pos, chromosome), snp_maf, cnv_pos, cnv_length, max_r2, p, distanceSNP_CNV(snp_pos, cnv_pos, cnv_length), distanceToTSS(snp_pos, chromosome)])
			snp_pos 	= int(line[1])
			cnv_pos 	= int(line[3])
			cnv_length 	= int(line[4])
			snp_maf 	= determineSNPMAF (float(line[5]))
			max_r2 		= float(line[7])
			p 		= float(line[8])

			hit_snp_pos.append(snp_pos)

	writeToOutputFile(output, chr)

	at_file = open('/data2/as/hitSNP-CNV/results/All.Del-20130807.LD-2000000.causal.pearson.chr' + str(chr) + '.txt')
	at_file.readline() # discard the header
	
	output = [] # clear output 
	
	line 		= at_file.readline().split()
	snp_pos 	= int(line[1])
	cnv_pos 	= int(line[3])
	cnv_length 	= int(line[4])
	snp_maf 	= determineSNPMAF (float(line[5]))
	max_r2 		= float(line[7])
	p 		= float(line[8])

	# loop over all association test results
	for test in at_file.readlines():
		line = [value for value in test.split()]
		if snp_pos == int(line[1]): # still same SNP?
			r2 = float(line[7])
			if not math.isnan(r2) and max_r2 < r2:
				max_r2 		= r2
				cnv_pos 	= int(line[3])
				cnv_length 	= int(line[4])
				p 		= float(line[8])
		else:
			if snp_pos not in hit_snp_pos: # not a hit-SNP?
				output.append([chr, 0, snp_pos, determineType(snp_pos, chromosome), snp_maf, cnv_pos, cnv_length, max_r2, p, distanceSNP_CNV(snp_pos, cnv_pos, cnv_length), distanceToTSS(snp_pos, chromosome)])
			
			snp_pos 	= int(line[1])
			cnv_pos 	= int(line[3])
			cnv_length 	= int(line[4])
			snp_maf 	= determineSNPMAF (float(line[5]))
			max_r2 		= float(line[7])
			p 		= float(line[8])

	writeToOutputFile(output, chr)

	


def main(): 

	for chr in range(1,23):
		chromosome = readChromosomeData(chr)
		processAT(chr, chromosome)

main()


