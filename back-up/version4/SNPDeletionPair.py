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
	Contains the functionality to process the intermediate results. 
"""

def readHitAndNonHitSNPDeletionPairs(results_file):
	hit_snp_deletion_pairs = [] 
	non_hit_snp_deletion_pairs = []
	for line in results_file:
		if len(line.split('\t')) < 18:
			continue 
		snp_deletion_pair = SNPDeletionPair(line)
		if snp_deletion_pair.snp.hit_snp:
			hit_snp_deletion_pairs.append(snp_deletion_pair)
		else:
			non_hit_snp_deletion_pairs.append(snp_deletion_pair)
	return hit_snp_deletion_pairs, non_hit_snp_deletion_pairs


def readSNPDeletionPairs(results_file, hit_snps_only = False, non_hit_snps_only = False):
	"""Returns all SNP-Deletion pairs in the given results_file."""
	snp_deletion_pairs = [] 

	if hit_snps_only:
		for line in results_file:
			if len(line.split('\t')) < 18:
				continue 
			snp_deletion_pair = SNPDeletionPair(line)
			if not snp_deletion_pair.hit_snp:
				continue
			snp_deletion_pairs.append(snp_deletion_pair)
	if non_hit_snps_only:
		for line in results_file:
			if len(line.split('\t')) < 18:
				continue 
			snp_deletion_pair = SNPDeletionPair(line)
			if snp_deletion_pair.hit_snp:
				continue
			snp_deletion_pairs.append(snp_deletion_pair)
	else:
		for line in results_file:
			snp_deletion_pairs.append(SNPDeletionPair(line))
	return snp_deletion_pairs

def strip(string):
	"""Strips the string on the left and right."""
	string = string.rstrip()
	return string.lstrip()


class SNP:
	"""Container for the data about the SNP."""
	def __init__(self, position, major_af, snp_type, dist_tss, hit_snp = False, hit_allele = None):
		self.position 	= position
		self.major_af 	= major_af
		self.type 	= snp_type
		self.dist_tss 	= dist_tss
		self.hit_snp  	= hit_snp
		self.hit_allele = hit_allele


class Deletion:
	"""Container for the data about the deletion."""
	def __init__(self, position, major_af, length):
		self.position 	= position
		self.major_af 	= major_af
		self.length	= length

class SNPDeletionPair: 
	"""Container for the data about a SNP-Deletion pair."""

	def __init__(self, line):
		"""Creates a results instance given the line in the file"""
		values 		= line.split('\t') # all values are stored here	
		self.chromosome	= values[0]
		self.hit_snp 	= False
		self.hit_allele = None
		if self.isHitSNP(values[1]): # determines whether this is a hit SNP or not
			self.hit_snp 	= True
			self.hit_allele = strip(values[2])
		
		# arguments: position, major_af, snp type ('intergenic'/'intronic'/'exonic'), distance to TSS
		self.snp = SNP(int(values[3]), self.readFloat(values[4]), strip(values[5]), int(values[6]), hit_snp = self.hit_snp, hit_allele = self.hit_allele)
		
		# arguments: position, major_af, length
		self.deletion = Deletion(int(values[7]), self.readFloat(values[10]), int(values[9])) # TODO change
		
		self.distance_snp_deletion = int(values[11])

		# 2x2 Table
		self.a = int(values[12])
		self.b = int(values[13])
		self.c = int(values[14])
		self.d = int(values[15])

		self.R = self.readFloat(values[16]) # Pearson's R
		self.p = self.readFloat(values[17]) # Fisher's p-value

	def print2(self): 
		print(self.chromosome, '\t', self.snp.position, '\t', self.snp.hit_allele, '\t', end = '') 
		print(self.deletion.position, '\t', self.deletion.length, '\t', end = '')
		print(self.R, '\t', self.p, '\t', self.a, '\t', self.b, '\t', self.c, '\t', self.d)

	def print(self):
		if self.hit_snp:
			print(self.chromosome, '\t', self.snp.position, '\t', self.snp.hit_allele, '\t', self.snp.type, '\t', self.snp.major_af, '\t', self.snp.dist_tss, '\t', self.distance_snp_deletion, '\t', end = '') # info on SNP
			print(self.deletion.position, '\t', self.deletion.length, '\t', self.deletion.major_af, '\t', end = '') # info on deletion
			print(self.a, '\t', self.b, '\t', self.c, '\t', self.d, '\t', self.R, '\t', self.p)
		
	
	def similarTo(self, other_snp_deletion_pair, major_af_threshold = 0.05, distance_snp_del_threshold = 10000, distance_tss_threshold = 10000):
		if self.snp.type != other_snp_deletion_pair.snp.type:
			return False
		if abs(self.snp.major_af - other_snp_deletion_pair.snp.major_af) > major_af_threshold:
			return False
		if abs(self.distance_snp_deletion - other_snp_deletion_pair.distance_snp_deletion) > distance_snp_del_threshold:
			return False 
		if abs(self.snp.dist_tss - other_snp_deletion_pair.snp.dist_tss) > distance_tss_threshold:
			return False 
		return True

	def materiallySignificant(self, r2_threshold = 0.8):
		if self.R**2.0 >= r2_threshold:
			return True
		return False	

	def isHitSNP(self, label):
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



