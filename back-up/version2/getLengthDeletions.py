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

import math 
import gzip
from scipy.stats.stats import pearsonr
import sys
from time import sleep

# Returns the major allele frequency (MAF) 
def determineMAF (af):
	if af >= .5: 
		return af
	else:
		return 1.0 - af


# Writes array to a prescripted output file 
def writeToOutputFile(output, chr):
	filename = 'length_deletions.txt' 
	output_file = open(filename, 'a')
	for length in output:
		output_file.write(str(chr) + '\t' + str(length)) 
		output_file.write('\n')

	output_file.close() 

# Returns the (zero-counting) index of the columns related to children 
def retrieveColumnsOfChildren (chr):
	# open the file with phased data
	phased_file = gzip.open('/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chr) + '_snps_and_dels.vcf.gz', 'r')

	# discard the headers (first 94 lines)
	for line in range(94):
		phased_file.readline() 

	columns_children = [] 

	header = [value for value in phased_file.readline().split()] 

	for i in range(len(header)):
		if header[i][-1] == 'c':
			columns_children.append(i)

	phased_file.close() 
	print 'Number of children chr' + str(chr) + ':',  len(columns_children)
	return columns_children


# Retrieves the deletions for a file:
# Returns:
#	deletions - every line a deletion
#		POS
#		LENGTH
#		Major allele frequency
#		Vector
def retrieveDeletions (chr, child):

	# open the file with phased data
	phased_file = gzip.open('/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chr) + '_snps_and_dels.vcf.gz', 'r')

	# discard the headers (first 94 lines)
	for line in range(95):
		phased_file.readline() 

	n_dels = 0 
	n_dels_with_maf1 = 0 

	deletions = [] # data is stored here

	# loop over all lines until you find a deletion
	for line in phased_file.readlines():
		l = [value for value in line.split()]
		if len(l[3]) > 1: # deletion! 
			# pos / length 
			deletions.append(len(l[3]) - 1)

	phased_file.close()
	print '---------------------------------' 
	print 'Chromosome ' + str(chr) + ': '
	print 'Deletions with major allele frequency < 1: ', n_dels
	print 'Deletions with major allele frequency = 1: ', n_dels_with_maf1
	print 'Total number of deletions: ', n_dels + n_dels_with_maf1
	print '---------------------------------' 
	return deletions



def main():
	for chr in range(3,22):
		child = retrieveColumnsOfChildren (chr)
		dels = retrieveDeletions(chr, child)
		writeToOutputFile(dels, chr)


main()
		
