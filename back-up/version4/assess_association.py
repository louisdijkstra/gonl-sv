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
import scipy as sp
import scipy.stats
from scipy.stats.stats import pearsonr

from DeletionReader import *
from SNPReader import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] <chromosome>

	<chromosome> - the chromosome that needs to be processed (1-22)

Processes the SNPs and deletions. For every SNP, it finds the deletions in a given range and determines Pearson's R and Fisher's p-value for association. 

Every line that is outputed represents a SNP-deletion pair. The columns represent:

	CHROM 		- the chromosome on which the SNP and deletion lie. 
	HIT		- is 'hit' when the SNP involved is a hit (GWAS) SNP, otherwise 'no_hit'. 
	HIT_ALLELE 	- if the SNP is a hit SNP, this is the hit allele. When the SNP is not a hit SNP, it is 'None'
	POS		- position of the SNP
	MAJOR_AF	- the major allele frequency of the SNP
	TYPE		- 'intergenic'/'intronic'/'exonic'
	DIST_TSS	- the distance to the nearest transcription start site (TSS)
	DEL_POS		- the position of the deletion as given in the VCF file
	DEL_START	- the start position of the deletion
	DEL_END		- the end of the deletion
	DEL_LENGTH	- the length of the deletion
	DEL_MAJOR_AF	- the major allele frequency of the deletion
	DIST_SNP_DEL	- the distance between the SNP and the deletion
	A		- 2x2 table entry
	B		- 2x2 table entry
	C		- 2x2 table entry
	D		- 2x2 table entry
	R		- Pearson's R
	p		- p-value resulting from applying Fisher's exact test to the 2x2 table 

The columns A,B,C and D provide the entire 2x2 table: 
			    Deletion	
		|	1	|	0	|
	--------|-------------------------------|-------
	      1	|	A	|	C	| A + C
	SNP   0	|	B	| 	D	| B + D
	--------|---------------------------------------
		|	A+B	|	C+D	| A+B+C+D
	
"""

class Table2x2:
	
	def __init__(self, classification1, classification2):
		self.a, self.b, self.c, self.d = 0, 0, 0, 0
		for i in range(len(classification1)):
			if classification1[i] == 1:
				if classification2[i] == 1: 	self.a += 1
				else:				self.b += 1
			else:
				if classification2[i] == 1:	self.c += 1
				else:				self.d += 1
		
	def R(self):
		"""Returns Pearson's R."""
		return float((self.a * self.d - self.b * self.c)) / math.sqrt((self.a+self.b)*(self.a+self.c)*(self.b+self.d)*(self.c+self.d))

	def FishersExactTest(self):
		[odds_ratio, p] = sp.stats.fisher_exact([[self.a,self.c],[self.b,self.d]],alternative='two-sided')
		return p 

	def store(self):
		print(self.a, '\t', self.b, '\t', self.c, '\t', self.d, '\t', end = '')

	

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("-a", action="store", dest="max_major_af", default=.96, type=float,
				  		help="Max. major allele frequency for a deletion and a snp to be considered. (Default = .96)")
	parser.add_option("-d", action="store", dest="distance_threshold", default=1000000, type=int,
				  		help="Max. distance between SNP and deletion for the pair to be considered. (Default = 10^6)")
	(options, args) = parser.parse_args()

	if (len(args)!=1):
		parser.print_help()
		return 1

	chromosome = int(args[0])

	deletion_reader = DeletionReader(max_major_af = options.max_major_af)
	# Read in the relevant data
	deletion_reader.read(chromosome)
	snp_reader = SNPReader(chromosome, max_major_af = options.max_major_af)		
	snp_reader.read()

	print('# of SNPs now: ', len(snp_reader.snps))

	# walk through all SNPs:
	for snp in snp_reader.snps:
		# walk through all deletions
		for deletion in deletion_reader.deletions: 
			distance_snp_deletion = deletion.minimumDistanceTo(snp.position)
			if distance_snp_deletion > options.distance_threshold:
				continue
			table = Table2x2(snp.haplotypes, deletion.haplotypes)
			# output relevant data:
			snp.store()
			deletion.store()
			print(distance_snp_deletion, '\t', end = '')
			table.store()
			print(table.R(), '\t', table.FishersExactTest())

if __name__ == '__main__':
	sys.exit(main())
