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

from HG19Reader import * 
from GWASCatalogueReader import * 
from HumanOmniReader import * 
from Variant import * # contains the functionality for dealing with SNPs and deletions

__author__ = "Louis Dijkstra"

"""
	Contains the reader for the phased SNPs.  
"""

class SNPReader: 
	"""Reads the deletions from file."""
	def __init__(self, chromosome, max_major_af = .96):
		self.max_major_af = max_major_af
	
		self.gwas_catalogue_reader = GWASCatalogueReader() 
		self.gwas_catalogue_reader.read(chromosome)

		self.human_omni_reader = HumanOmniReader()
		self.human_omni_reader.read(chromosome) 

		self.hg19_reader = HG19Reader()
		self.hg19_reader.read(chromosome)

		self.vcf_reader = vcf.Reader(open(self.returnLocation(chromosome)))

	def returnLocation(self, chromosome):
		if chromosome != 1 and chromosome != 22:
			return '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chromosome) + '_snps_and_dels.vcf.gz'
		elif chromosome == 1:
			return '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr1_snps_and_dels.vcf'
		return '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr22_snps_and_dels.vcf'
	
	def readHitSNPs(self):
		"""Reads in only the hit SNPs present in the given VCF file."""
		self.snps = []
		for vcf_record in self.vcf_reader: 
			if not isSNP(vcf_record): # if vcf is not a SNP
				continue
			snp = SNP(vcf_record)	
			if snp.major_af > max_major_af: # major allele frequency too high
				continue
			if not self.human_omni_reader.snpPresent(snp.position):  # not present in the HumanOmni data set
				continue
			hit_allele = self.gwas_catalogue_reader.getHitAllele(snp.position)
			if hit_allele == None:
				continue
			
			snp_type, dist_tss = self.hg19_reader.determineTypeAndDistanceTSS(snp.position)
			self.snps.append(SNP(vcf_record, snp_type = snp_type, hit_allele = hit_allele, dist_tss = dist_tss))			


	def readRawList(self, snp_positions): 
		"""Reads in all SNPs for which their position is present in the given list.
		   Their distance to TSS, type and (when relevant) hit allele are not determined."""
		self.snps = []
		for vcf_record in self.vcf_reader:
			if not isSNP(vcf_record):
				continue
			if not vcf_record.POS in snp_positions:
				continue
			snp = SNP(vcf_record)	
			if snp.major_af > self.max_major_af:
				continue
			self.snps.append(snp)
			
			
	def read(self):
		"""Reads all the SNPs present in the given VCF file. 
 		   NOTE: this function must be called first before any other function can be called."""
		self.snps = [] 
		for vcf_record in self.vcf_reader: 
			if not isSNP(vcf_record):
				continue
			snp = SNP(vcf_record)	
			if snp.major_af > self.max_major_af:
				continue
			if not self.human_omni_reader.snpPresent(snp.position): 
				continue

			snp_type, dist_tss = self.hg19_reader.determineTypeAndDistanceTSS(snp.position)
			hit_allele = self.gwas_catalogue_reader.getHitAllele(snp.position)
			self.snps.append(SNP(vcf_record, snp_type = snp_type, hit_allele = hit_allele, dist_tss = dist_tss))			
		
	def getHitSNPs(self):
		self.hit_snps = []
		for snp in self.snps:
			if snp.isHitSNP():
				self.hit_snps.append(snp)
			

	def printAll(self):
		for snp in snps:
			snp.print()
		
 
