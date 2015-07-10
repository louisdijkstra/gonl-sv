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

usage = """%prog [options] <snp-deletion-pair-file1> ... <snp-deletion-pair-file_m1>

Reads in multiple '.pairs' files and returns the 

- percentage of items deemed significant 
- a list of significant pairs 
"""

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("--children", action="store_true", dest="children_only", default=False, 
					help = "Extract the haplotypes of the children only. In case of twins, only the first is selected.")
	parser.add_option("--indels", action="store_true", dest="indels_only", default=False, 
				help = "Extract the haplotypes of the indels only")
	parser.add_option("--parents", action="store_true", dest="parents_only", default=False, 
					help = "Extract the haplotypes of the parents only")
	parser.add_option("--snps", action="store_true", dest="snps_only", default=False, 
				help = "Extract the haplotypes of the SNPs only")
	parser.add_option("-a", action="store", dest="minor_af_thres", default=0.0, type=float,
				help = "Extract the haplotypes of those variants for which the minor allele frequency is equal or larger than this value. (Default = 0.0)")
	(options, args) = parser.parse_args()

	if (len(args)!=1):
		parser.print_help()
		return 1

	vcf_reader = vcf.Reader(open(args[0]))
	
	# allocate memory
	autosome, ones, n = 0, 0, 0 
	type_variant = "*"
	haplotype = []
	allele_freq = 0.0
	if options.parents_only:
		for vcf_record in vcf_reader: 
			# determine type 
			type_variant = '*'
			if 'SVTYPE' in vcf_record.INFO:
				if isinstance(vcf_record.INFO['SVTYPE'], list):
					if vcf_record.INFO['SVTYPE'][0] == 'DEL':
						type_variant = '-'
					elif vcf_record.INFO['SVTYPE'][0] == 'INS':
						type_variant = '+'
			else:
				if len(vcf_record.REF) > 1 and len(vcf_record.ALT[0]) == 1:
					type_variant = '-'
				elif len(vcf_record.REF) == 1 and len(vcf_record.ALT[0]) > 1:
					type_variant = '+'

			# check whether type is under consideration
			if type_variant == '*' and options.indels_only:
				continue
			elif type_variant != '*' and options.snps_only:
				continue


			n, ones = 0, 0
			haplotype = []
			for s in vcf_record.samples:
				if s.sample[-1] == 'a' or s.sample[-1] == 'A' or s.sample[-1] == 'b' or s.sample[-1] == 'B':
					haplotype.append(int(s['GT'][0]))
					haplotype.append(int(s['GT'][2]))
					ones += int(s['GT'][0]) + int(s['GT'][2])
					n += 2 
			allele_frequency = float(ones)/float(n)
			
			if allele_frequency < options.minor_af_thres or 1.0 - allele_frequency < options.minor_af_thres:
				continue

			if type_variant == '*' :
				print("%d %c %d . "%(returnAutosome(vcf_record.CHROM), type_variant, vcf_record.POS), end='')
				print(vcf_record.REF, vcf_record.ALT[0], end=' ')
			elif type_variant == '-':
				print("%d %c %d %d . . "%(returnAutosome(vcf_record.CHROM), type_variant, vcf_record.POS, returnIndelLength(vcf_record)), end='')
			elif type_variant == '+':
				print("%d %c %d %d . . "%(returnAutosome(vcf_record.CHROM), type_variant, vcf_record.POS, returnIndelLength(vcf_record)), end='')

			print("%d %f "%(n, allele_frequency), end='')
			for i in range(n):
				print('%d'%haplotype[i], end='')
			print('')	
	elif options.children_only:
		for vcf_record in vcf_reader: 
			# determine type 
			type_variant = '*'
			if 'SVTYPE' in vcf_record.INFO:
				if isinstance(vcf_record.INFO['SVTYPE'], list):
					if vcf_record.INFO['SVTYPE'][0] == 'DEL':
						type_variant = '-'
					elif vcf_record.INFO['SVTYPE'][0] == 'INS':
						type_variant = '+'
			else:
				if len(vcf_record.REF) > 1 and len(vcf_record.ALT[0]) == 1:
					type_variant = '-'
				elif len(vcf_record.REF) == 1 and len(vcf_record.ALT[0]) > 1:
					type_variant = '+'

			# check whether type is under consideration
			if type_variant == '*' and options.indels_only:
				continue
			elif type_variant != '*' and options.snps_only:
				continue


			n, ones = 0, 0
			haplotype = []
			for s in vcf_record.samples:
				if s.sample[-1] == 'c' or s.sample[-1] == 'C':
					haplotype.append(int(s['GT'][0]))
					haplotype.append(int(s['GT'][2]))
					ones += int(s['GT'][0]) + int(s['GT'][2])
					n += 2 
			allele_frequency = float(ones)/float(n)
			
			if allele_frequency < options.minor_af_thres or 1.0 - allele_frequency < options.minor_af_thres:
				continue

			if type_variant == '*' :
				print("%d %c %d . "%(returnAutosome(vcf_record.CHROM), type_variant, vcf_record.POS), end='')
				print(vcf_record.REF, vcf_record.ALT[0], end=' ')
			elif type_variant == '-':
				print("%d %c %d %d . . "%(returnAutosome(vcf_record.CHROM), type_variant, vcf_record.POS, returnIndelLength(vcf_record)), end='')
			elif type_variant == '+':
				print("%d %c %d %d . . "%(returnAutosome(vcf_record.CHROM), type_variant, vcf_record.POS, returnIndelLength(vcf_record)), end='')

			print("%d %f "%(n, allele_frequency), end='')
			for i in range(n):
				print('%d'%haplotype[i], end='')
			print('')
	else:
		for vcf_record in vcf_reader: 
			# determine type 
			type_variant = '*'
			if 'SVTYPE' in vcf_record.INFO:
				if isinstance(vcf_record.INFO['SVTYPE'], list):
					if vcf_record.INFO['SVTYPE'][0] == 'DEL':
						type_variant = '-'
					elif vcf_record.INFO['SVTYPE'][0] == 'INS':
						type_variant = '+'
			else:
				if len(vcf_record.REF) > 1 and len(vcf_record.ALT[0]) == 1:
					type_variant = '-'
				elif len(vcf_record.REF) == 1 and len(vcf_record.ALT[0]) > 1:
					type_variant = '+'

			# check whether type is under consideration
			if type_variant == '*' and options.indels_only:
				continue
			elif type_variant != '*' and options.snps_only:
				continue


			n, ones = 0, 0
			haplotype = []
			for s in vcf_record.samples:
				haplotype.append(int(s['GT'][0]))
				haplotype.append(int(s['GT'][2]))
				ones += int(s['GT'][0]) + int(s['GT'][2])
				n += 2 
			allele_frequency = float(ones)/float(n)
			
			if allele_frequency < options.minor_af_thres or 1.0 - allele_frequency < options.minor_af_thres:
				continue

			if type_variant == '*' :
				print("%d %c %d . "%(returnAutosome(vcf_record.CHROM), type_variant, vcf_record.POS), end='')
				print(vcf_record.REF, vcf_record.ALT[0], end=' ')
			elif type_variant == '-':
				print("%d %c %d %d . . "%(returnAutosome(vcf_record.CHROM), type_variant, vcf_record.POS, returnIndelLength(vcf_record)), end='')
			elif type_variant == '+':
				print("%d %c %d %d . . "%(returnAutosome(vcf_record.CHROM), type_variant, vcf_record.POS, returnIndelLength(vcf_record)), end='')

			print("%d %f "%(n, allele_frequency), end='')
			for i in range(n):
				print('%d'%haplotype[i], end='')
			print('')

if __name__ == '__main__':
	sys.exit(main())

