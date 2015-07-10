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
	Contains the classes Variation and the subclases SNP and Deletion. Used for processing VCF files. 
"""


# TODO!

#def snpsSimilar (snp1, snp2, dist_tss_threshold = 1000, dist_snp_del_threshold = 1000):
#	if snp1.type != snp2.type:
#		return False
#	if abs(snp1.dist_tss - snp2.dist_tss) > dist_tss_threshold:
#		return False
	

def isDeletion (vcf_record):
	"""Determines whether a vcf record represents a deletion or not. NOTE: only considers first alternative (ALT[0]); others are neglected"""
	if 'SVTYPE' in vcf_record.INFO: # checks whether the key SVTYPE is present in INFO
		if vcf_record.INFO['SVTYPE'][0] == 'DEL':	
			return True
	else: # makes statement on the basis of length of ALT and REF 
		if len(vcf_record.REF) > 1 and len(vcf_record.ALT[0]) == 1:
			return True
	return False

def isSNP (vcf_record):
	"""Determines whether a vcf record represents a SNP or not. NOTE: only considers first alternative (ALT[0]); others are neglected"""
	if len(vcf_record.REF) == 1 and len(vcf_record.ALT[0]) == 1:
		return True
	return False


def relevantSample(sample):
	"""Returns true when the sample is relevant, and false otherwise. If it ends, for example, on a 'c', than it is a kid etc."""
	ending = sample[-1]

	if ending == 'a' or ending == 'A' or ending == 'b' or ending == 'B': # only consider parents
		return True
	return False


class SNPDeletionPair:
	
	def __init__(self, snp, deletion):
		self.snp = snp
		self.deletion = deletion 
		

class Variant: 
	"""Class for representing a variation given a VCF record."""
	def __init__(self, vcf_record):
		self.vcf_record = vcf_record 
		self.chromosome = str(vcf_record.CHROM)
		self.position 	= vcf_record.POS 
		self.determineHaplotypes()
		self.determineMajorAlleleFrequency()	

	def determineHaplotypes(self):
		self.haplotypes = []
		for sample in self.vcf_record.samples:		
			if not relevantSample(sample.sample):
				continue
			chromosome1 = int(sample['GT'][0])
			chromosome2 = int(sample['GT'][2])
			self.haplotypes.append(chromosome1)
			self.haplotypes.append(chromosome2)

	def determineMajorAlleleFrequency(self):
		"""Determines the major allele frequency on the basis of the haplotypes."""
		self.major_af = sum(self.haplotypes) / float(len(self.haplotypes))
		if self.major_af < 0.5:
			self.major_af = 1.0 - self.major_af
		
class Deletion(Variant):
	"""Class for representing a deletion."""
	def __init__(self, vcf_record):
		Variant.__init__(self, vcf_record)
		self.start = vcf_record.POS + 1
		self.length = abs(len(vcf_record.REF) + 1 - len(vcf_record.ALT[0]))
		self.end   = self.start + self.length - 1

	def print(self):
		"""Prints info on the deletion to the command line."""
		print('Deletion\t', self.chromosome, '\t', self.position, '\t', self.length, '\t', self.haplotypes)

	def minimumDistanceTo(self, location):
		"""Returns the minimum distance to a given location."""
		if location < self.start:
			return self.start - location
		elif location > self.end:
			return location - self.end
		return 0 

	def store(self):
		"""Used for printing info for storing."""
		print(self.position, '\t', self.start, '\t', self.end, '\t', self.length, '\t', self.major_af, '\t', end = '')

class SNP(Variant):
	"""Class for representing a SNP."""
	def __init__(self, vcf_record, snp_type = None, hit_allele = None, dist_tss = None):
		Variant.__init__(self, vcf_record)	
		self.type 	= snp_type
		self.hit_allele = hit_allele
		self.dist_tss 	= dist_tss # distance to the nearest transcription start site (TSS)
		if hit_allele == None:
			self.hit_snp = False
		else:
			self.hit_snp = True
			self.checkHaplotypes()

	def checkHaplotypes(self):
		"""When the SNP is a hit SNP, it might be that the haplotypes must be flipped, from 0 to 1 and the other way around."""
		if self.vcf_record.REF != self.hit_allele:
			for i in range(len(self.haplotypes)):
				if self.haplotypes[i] == 0:
					self.haplotypes[i] = 1
				else:
					self.haplotypes[i] = 0

	def isHitSNP(self):
		return self.hit_snp

	def sameType (self, snp_type):
		return self.type == snp_type

	def print(self):
		"""Prints info on the SNP to the command line."""
		print('SNP\t\t', self.chromosome, '\t', self.position, '\t', self.type, '\t', self.hit_allele, '\t', self.dist_tss, '\t', self.haplotypes)


	def store(self):
		"""Used for printing info for storing."""
		print(self.chromosome, '\t', end = '')
		if self.hit_snp:
			print('hit\t', self.hit_allele, '\t', end = '')
		else:
			print('no_hit\tNone\t', end = '')
		print(self.position, '\t', self.major_af, '\t', self.type, '\t', self.dist_tss, end = '')



