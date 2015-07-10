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
import random
import scipy as sp
import scipy.stats
from scipy.stats.stats import pearsonr

from SNPReader import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] <tag_snp_file>

Adds the reference and alternative allele to the tag SNP data. 

The output is presented as a table. Every row is for one deletion. The columns (seperated by a tab) represent: 

	CHR 		- autosome (integer)
	POS		- position of the deletion 
	LENGTH		- length of the deletion 
	SNP_POS		- position of the tag SNP 
	ALLELE #1	- the allele of the tag SNP with the higher allele frequency
	ALLELE #2	- the allele of the tag SNP with the lower allele frequency
	A		- entry of the 2x2 contingency table (see below)
	B		- entry of the 2x2 contingency table (see below)
	C		- entry of the 2x2 contingency table (see below)
	D		- entry of the 2x2 contingency table (see below)
	R		- Pearson's R computed on the basis of the contingency table
	p 		- p-value of Fisher's exact test the basis of the contingency table
	P(DEL|ALLELE #1)- the estimated probability of observing a deletion given the presence of allele #1
	P(DEL|ALLELE #2)- the estimated probability of observing a deletion given the presence of allele #2

The columns A, B, C and D provide all information about the 2x2 contigency table: 

				     	             Deletion
				|	Present		|	Absent 		|	total
	----------------------------------------------------------------------------------------
	allele #1		|	   A		|	  C		|     	A + C
SNP	allele #2		|	   B		|	  D		|	B + D
	----------------------------------------------------------------------------------------
	total			| 	  A + B		| 	 C + D		|      A+B+C+D
"""

class DeletionTAGSNPPair:
	
	def __init__(self, line):
		self.line = line # store original
		values = line.split('\t') # all values are stored here
		self.chromosome = int(values[0])
		self.pos 	= int(values[1])
		self.length	= int(values[2])
		self.snp_pos 	= int(values[3])
		self.a		= int(values[4])
		self.b		= int(values[5])
		self.c		= int(values[6])
		self.d		= int(values[7])
		self.R		= float(values[8])
		self.p		= float(values[9])

	def addReferenceAlternativeAllele(self, ref, alt):
		self.ref = ref
		self.alt = alt

	def update(self):
		if self.b + self.d > self.a + self.c: # turn reference and alternative allele around
			self.R = -1.0 * self.R
			temp = self.ref
			self.ref = self.alt
			self.alt = temp 
			temp = self.b
			self.b = self.a
			self.a = temp
			temp = self.d
			self.d = self.c
			self.c = temp
			
	def print(self):
		# print data related to the deletion:
		print(self.chromosome, '\t', self.pos, '\t', self.length, '\t', end = '')
		# print data related to the SNP:
		print(self.snp_pos, '\t', self.ref, '\t', self.alt, '\t', end = '')		
		# print data about the table 
		print(self.R, '\t', self.p, '\t', self.a, '\t', self.b, '\t', self.c, '\t', self.d, '\t', end = '')
		# compute conditional probabilities and print them
		print(float(self.a) / float(self.a + self.c), '\t', float(self.b) / float(self.b + self.d))


def main():
	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()

	list_of_autosomes = range(1,23)

	for autosome in list_of_autosomes: 
		
		# read tag SNP file 
		tag_snp_file = open("Results/tags_chr" + str(autosome) + ".txt")
		tag_snp_file.next() # discard header

		#print('Processing tag SNP file...')
		del_snp_pairs, snp_positions = [], []
		for line in tag_snp_file:
			del_snp_pairs.append(DeletionTAGSNPPair(line))
			snp_positions.append(del_snp_pairs[-1].snp_pos)
		#print('DONE processing tag SNP file...')
	
		# read VCF file
		#print('Processing VCF file...')
		snp_reader = SNPReader(1, max_major_af = .96)		
		snp_reader.readRawList(snp_positions)
		#print('DONE Processing VCF file...')
		
		# Postprocess every tag SNP-deletion pair
		for del_snp_pair in del_snp_pairs:
			for snp in snp_reader.snps:
				if del_snp_pair.snp_pos == snp.position: 
					del_snp_pair.addReferenceAlternativeAllele(snp.vcf_record.REF, snp.vcf_record.ALT[0])
			del_snp_pair.update()
			del_snp_pair.print()	

if __name__ == '__main__':
	sys.exit(main())

