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

from Variant import * # contains the functionality for dealing with deletions

__author__ = "Louis Dijkstra"

"""
	Contains the reader for the phased deletions. 
"""

class DeletionReader: 
	"""Reads the deletions from file."""
	def __init__(self, location = "/data1/gonl/new-deletion-phasing/deletions-20-10000-phased-20140925.vcf", max_major_af = .96):
		self.vcf_reader = vcf.Reader(open(location))
		self.max_major_af = max_major_af
		self.deletions = [] 

	def read(self, chromosome):
		"""Reads all the deletions on one chromosome. This must be an integer. 
 		   NOTE: this function must be called first before any other function can be called."""
		self.deletions = [] # reset
		chromosome = str(chromosome)
		#n_deletions_processed = 0 
		for vcf_record in self.vcf_reader: 
			if vcf_record.CHROM != chromosome:
				continue
			#n_deletions_processed += 1
			#if n_deletions_processed % 1000 == 0:
			#	print(n_deletions_processed, ' deletions processed so far.')
			deletion = Deletion(vcf_record)
			if deletion.major_af > self.max_major_af:
				continue
			self.deletions.append(Deletion(vcf_record))
		
	def printAll(self):
		for deletion in deletions:
			deletion.print()
	
		
 
