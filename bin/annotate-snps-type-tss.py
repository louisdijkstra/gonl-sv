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

from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import math
import vcf

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))[:-3] + 'python')
from Variant import *

__author__ = "Louis Dijkstra"

usage = """%prog <haplotype-file> <hg19-file>

Reads in a haplotype file <haplotype-file> as generated by 
'extract-haplotypes.py' and annotates every SNP in the file 
with its type (intergenic/intronic/exonic) and the distance
to the nearest TSS. 
"""

def main():
	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()

	if (len(args)!=2):
		parser.print_help()
		return 1

	haplo_file = open(args[0], 'r')
	hg19_file  = open(args[1], 'r')

	# data from hg19 is stored here
	hg19_n_genes	= [0 for i in range(1,23)]
	hg19_tss  	= [[] for i in range(1,23)] # transcription start site
	hg19_scr  	= [[] for i in range(1,23)] # start coding region
	hg19_ecr 	= [[] for i in range(1,23)] # end coding region
	hg19_n_exons  	= [[] for i in range(1,23)] 
	hg19_exon_start = [[] for i in range(1,23)]
	hg19_exon_end 	= [[] for i in range(1,23)]

	# skip header
	hg19_file.next()
	
	# read in and process the hg19 file
	for line in hg19_file:
		values = line.split()
		#print(values)
		index = returnAutosome(values[2]) - 1 
		if index >= 22: # sex chromosome or worse! 
			continue 
		hg19_n_genes[index] 	+= 1 
		hg19_tss[index].append(int(values[4]))
		hg19_scr[index].append(int(values[6]))
		hg19_ecr[index].append(int(values[7]))

		# read in the exon starts and ends
		hg19_n_exons[index].append(int(values[8]))
		exon_starts, exon_ends = [], []
		for exon_start in values[9][:-1].split(','):
			exon_starts.append(int(exon_start)) 
		for exon_end in values[10][:-1].split(','):
			exon_ends.append(int(exon_end)) 

		hg19_exon_start[index].append(exon_starts)
		hg19_exon_end[index].append(exon_ends)
		
	hg19_file.close()

	# allocate memory
	values, autosome, type_variant, pos, length, ref, alt, n, allele_freq, haplotypes, allele, dist_tss, snp_type = [], 0, '*', 0, 0, '.', '.', 0, 0.0, [], '.', 0, '.'

	for line in haplo_file:
		# read in the line
		values 		= line.split()
		autosome 	= int(values[0]) 
		index 		= autosome - 1
		type_variant 	= values[1].strip()
		pos 		= int(values[2])
		length 		= values[3].strip()
		ref 		= values[4].strip()
		alt 		= values[5].strip()
		n 		= int(values[6])
		allele_freq 	= float(values[7])
		haplotypes 	= values[8].strip()
		
		if len(values) < 10:
			allele = '?'
		else:
			allele = values[9].strip()

		if type_variant == '+' or type_variant == '-': # in case of an indel
			print('%d %c %d '%(autosome, type_variant, pos), end='')
			print(length, ref, alt, n, allele_freq, haplotypes, allele, '.', '.') 
			continue 
		
		# determine distance to nearest TSS
		dist_tss = abs(pos - hg19_tss[index][0])
		for i in range(1, len(hg19_tss[index])):
			if abs(pos - hg19_tss[index][i]) < dist_tss:
				dist_tss = abs(pos - hg19_tss[index][i]) 

		snp_type = 'intergenic'
		
		for i in range(hg19_n_genes[index]):
			if hg19_scr[index][i] <= pos <= hg19_ecr[index][i]: # in coding region
				# determine whether it is in an exon
				for j in range(hg19_n_exons[index][i]):
					if hg19_exon_start[index][i][j] <= pos <= hg19_exon_end[index][i][j]:
						snp_type = 'exonic' 
						break 
					
				snp_type = 'intronic'
				
				break 

		print('%d %c %d '%(autosome, type_variant, pos), end='')
		print(length, ref, alt, n, allele_freq, haplotypes, allele, dist_tss, snp_type) 

		
	
if __name__ == '__main__':
	sys.exit(main())

