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

__author__ = "Louis Dijkstra"

"""
	The class HG19Reader helps with reading the hg19 file. It aids in determining whether a SNP is 
	'intronic', 'exonic' or 'intergenic'. 
"""

class Exon:
	"""Simple class to represent an exon."""
	def __init__(self, start, end):
		self.start = start
		self.end   = end

	def contains(self, location):
		"""Returns true when the location lies inside the exon, otherwise false."""
		if location >= self.start and location <= self.end:
			return True
		return False

class Gene: 
	"""Class to represent the data coming from the hg19 file."""
	def __init__(self, raw_data):
		self.chromosome 		= raw_data[2]
		self.transcription_start_site 	= int(raw_data[4])
		self.transcription_end_site 	= int(raw_data[5])
		self.start_coding_region 	= int(raw_data[6])
		self.end_coding_region 		= int(raw_data[7])
		self.n_exons 			= int(raw_data[8]) # number of exons
		# read exons		
		self.exons 	= []
		exon_starts	= raw_data[9].split(',')
		exon_ends 	= raw_data[10].split(',')
		for i in range(self.n_exons):
			self.exons.append(Exon(exon_starts[i], exon_ends[i]))

	def returnDistanceToTSS(self, location):
		"""Returns the distance to the transcription start site (TSS)."""
		return abs(self.transcription_start_site - location)

	def print(self):
		print(self.chromosome, '\t', self.transcription_start_site, '\t', self.transcription_end_site, '\t', self.start_coding_region, '\t', self.end_coding_region, self.n_exons)

class HG19Reader:

	def __init__(self, location = '/ufs/dijkstra/Projects/SNPs_CNVs/data/hg19.txt'):
		self.file = open('/ufs/dijkstra/Projects/SNPs_CNVs/data/hg19.txt', 'r') 
		self.genes = [] # list of genes. List is filled when the file is read

	def read(self, chromosome):
		"""Reads all genes on one chromosome. This must be an integer. 
 		   NOTE: this function must be called first before any other function can be called."""
		self.file.readline() # header must be discarded
		self.genes = [] # reset 
		chromosome = 'chr' + str(chromosome)
		for line in self.file:
			data = [value for value in line.split()]
			if data[2] != chromosome:
				continue
			self.genes.append(Gene(data)) 


	def determineTypeAndDistanceTSS(self, snp_position):
		"""Returns both the type and the distance to the nearest transcription start site (TSS)."""
		min_dist = float('inf')
		snp_type = 'intergenic'
		for gene in self.genes: 
			distance_to_tss = gene.returnDistanceToTSS(snp_position)
			if min_dist > distance_to_tss:
				min_dist = distance_to_tss
			# type 
			if gene.start_coding_region <= snp_position <= gene.end_coding_region:
				# determine whether location is in an exon
				for exon in gene.exons:
					if exon.contains(snp_position):
						snp_type = 'exonic'
				snp_type = 'intronic'
 		return snp_type, min_dist

	def returnMinimumDistanceToTSS(self, location):
		"""Returns the minimal distance to a transcription start site (TSS)."""
		min_dist = float('inf')
		for gene in self.genes:
			distance_to_tss = gene.returnDistanceToTSS(location)
			if min_dist < distance_to_tss:
				min_dist = distance_to_tss
		return min_dist

	def determineType(self, location):
		"""Returns the type of a particular location: either 'intergenic', 'intronic' or 'exonic'."""
		for gene in self.genes:
			if gene.start_coding_region <= location <= gene.end_coding_region:
				# determine whether location is in an exon
				for exon in gene.exons:
					if exon.contains(location):
						return 'exonic'
				return 'intronic'
		return 'intergenic'

	def print(self):
		for gene in self.genes:
			gene.print()

	def close(self):
		self.file.close()
