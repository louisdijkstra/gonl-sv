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
	The class GWASCatalogueReader helps with reading the GWAS catalogue file. 
"""

class HitSNP:
	"""Class for representing Hit SNPs from the GWAS catalogue."""
	def __init__(self, position, hit_allele):
		self.position = position
		self.hit_allele = hit_allele

	def print(self):
		print(self.position, '\t', self.hit_allele)

class GWASCatalogueReader:
	
	def __init__(self, location = '/ufs/dijkstra/Projects/SNPs_CNVs/data/GWAS_catalogue_caucasian_SNPIDs.alleles.txt'):
		self.file = open(location, 'r') 
		self.hit_snps = [] # list of hit snps. List is filled when the file is read
	
	def read(self, chromosome):
		"""Reads all the hit SNPs on one chromosome. This must be an integer. 
 		   NOTE: this function must be called first before any other function can be called."""
		self.hit_snps = [] # reset
		chromosome = 'chr' + str(chromosome)
		for line in self.file.readlines():
			raw_data = [value for value in line.split()]
			if raw_data[0] != chromosome:
				continue
			if len(raw_data) != 4: # data is not complete 
				continue 
			self.hit_snps.append(HitSNP(int(raw_data[1]), raw_data[3]))

	def getHitAllele(self, position):
		"""Returns the hit allele from a location. If not present in the list, returns 'None'."""
		for hit_snp in self.hit_snps:
			if hit_snp.position != position:
				continue
			return hit_snp.hit_allele
		return None

	def print(self):
		for hit_snp in self.hit_snps:
			hit_snp.print()
			
	def close(self):
		self.file.close()
