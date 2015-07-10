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

# hg19.py is used to read chromosomal data


# Returns the type of a location (intergenic/intronic/exonic)given 
# the chromosomal data (chromosome)
def determineType (position, genes):
	# move through all coding regions in the chromosome
	for gene in genes:
		if gene['start_coding_region'] <= position <= gene['end_coding_region']:
			# determine whether SNP is in an exon
			for exon in range(gene['n_exons']):
				if gene['exons'][exon][0] <= position <= gene['exons'][exon][1]:
					return 2
			return 1
	return 0

# Returns the distance to the nearest TSS given chromosomal data (chromosome)
def distanceToTSS (position, genes):
	min_dist = float('inf') 
	for gene in genes:
		dist = abs(gene['tss'] - position)
		if dist < min_dist:
			min_dist = dist 
	return min_dist 

def print_hg19 (genes):
	print 'tss\ttes\tstart_coding_region\tend_coding_region\tn_exons\texons' 
	for gene in genes:
		print gene['tss'], '\t', gene['tes'], '\t', gene['start_coding_region'], '\t', gene['end_coding_region'], '\t', gene['n_exons'], '\t', gene['exons']

# Reads chromosomal data from file for a specified chromosome
def read (chromosome):	
	hg19_file = open('/ufs/dijkstra/Projects/SNPs_CNVs/data/hg19.txt', 'r') # open data file
	hg19_file.readline() # discard the header 

	genes = [] # data is stored here

	# Read through all lines...
	for line in hg19_file.readlines():
		data = [value for value in line.split()]
		if data[2] == 'chr' + str(chromosome):
			gene = [('tss', int(data[4])),
				('tes', int(data[5])),
				('start_coding_region', int(data[6])),
				('end_coding_region', int(data[7])),
				('n_exons', int(data[8]))] 
			exon_starts 	= data[9].split(',')
			exon_ends	= data[10].split(',')
			gene.append(('exons', zip(map(int, exon_starts[:int(data[8])]), map(int, exon_ends[:int(data[8])]))))			
			genes.append(dict(gene))
	hg19_file.close()
	return genes 


