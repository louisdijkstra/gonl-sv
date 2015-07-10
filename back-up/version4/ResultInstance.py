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

__author__ = "Louis Dijkstra"

"""
	Contains the functionality to process the intermediate results. 
"""

def strip(string):
	"""Strips the string on the left and right."""
	string = string.rstrip()
	return string.lstrip()


class SNP:
	"""Container for the data about the SNP."""
	def __init__(self, position, major_af, snp_type, dist_tss, hit_snp = False, hit_allele = None):
		self.position 	= position
		self.major_af 	= major_af
		self.snp_type 	= snp_type
		self.dist_tss 	= dist_tss
		self.hit_snp  	= hit_snp
		self.hit_allele = hit_allele

class Deletion:
	"""Container for the data about the deletion."""
	def __init__(self, position, major_af, start, end, length):
		self.position 	= position
		self.major_af 	= major_af
		self.start 	= start
		self.end 	= end
		self.length	= lenght

class ResultInstance: 
	"""Container for the data about a SNP-Deletion pair."""

	def __init__(self, line):
		"""Creates a results instance given the line in the file"""
		self.values 		= line.split('\t') # all values are stored here	
		self.chromosome		= self.values[0]
		if self.isHitSNP(self.values[1]) # determines whether this is a hit SNP or not
			self.hit_snp 	= True
			self.hit_allele = strip(self.values[2])
		else:
			self.hit_snp 	= False
			self.hit_allele = None
		self.position 		= int(self.values[3])
		self.major_af		= 
				
		
		

	def isHitSNP(label):
		if strip(label) == 'hit':
			return True
		return False

	def readFloat(self, x):
		x = x.rstrip()
		x = x.lstrip()
		if x == "None":
			return None
		elif x == "-inf":
			return float("-inf")
		else:
			return float(x)



