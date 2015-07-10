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

# vcf.py contains various functions for reading and (partly) processing
#	.vcf file formats

import gzip 		# helps to read .gz 

# open_vcf_file opens a vcf file and returns the file pointer
def openVCFFile (vcf_filename): 
	vcf_file = None # allocate file pointer
	if '.gz' in vcf_filename:
		vcf_file = gzip.open(vcf_filename, 'r')
	else:
		vcf_file = open(vcf_filename, 'r')
	return vcf_file 

# discards the header	# TODO check the number of lines that need discarding 
def discardVCFHeaders (vcf_file, number_of_lines=93):
	for line in range(number_of_lines+1): vcf_file.readline() 

# returnColumnsOfChildren returns the indices of the columns 
# related to children. 
def returnColumnsOfChildren (vcf_filename):
	vcf_file = openVCFFile(vcf_filename)
	discardVCFHeaders(vcf_file)
	indices_children = [] 
	index = 0 
	for column_name in vcf_file.readline().split():
		if column_name[-1] == 'c':
			indices_children.append(index)
		index = index + 1 
	vcf_file.close()
	return indices_children 


def returnPhaseWithAllele(line, indices, hitallele):
	phase = returnPhase(line, indices)
	if line[4] == hitallele:
		return phase
	elif line[3] == hitallele:
		for i in range(len(phase)):
			phase[i] = 1 - phase[i]
		return phase
	else:
		print 'WAIT!? the hit allele does not appear for this SNP...'
		return phase 

def returnPhase(line, indices):
	phase = [] 
	for index in indices:
		phase.append(int(line[index].split('|')[0]))
		phase.append(int(line[index].split('|')[1]))
	return phase


def determineAlleleFrequency (phase):
	af = sum(phase) / float(len(phase))
	return af
 
def determineMajorAlleleFrequency (phase):
	af = sum(phase) / float(len(phase)) 
	if af < .5:
		return 1 - af
	else:
		return af 	

def printDeletions (deletions):
	print 'Position\tLength\tMajor AF\tPhase'
	for d in deletions:
		print d['pos'], '\t', d['length'], '\t', d['af'], '\t', d['phase']

def returnDeletions(chr, indices, majorAF_threshold = 1):
	vcf_file = openVCFFile('/data2/as/hitSNP-CNV/data/phased-deletions/Chr' + str(chr) + '.vcf')
	discardVCFHeaders(vcf_file, number_of_lines=5) 
	
	deletions = [] # data is stored here
	
	for line in vcf_file.readlines():
		info = [value for value in line.split()]
		deletion = [('pos', int(info[1])), ('length', len(info[3]) - 1)]
		phase = returnPhase(info, indices)
		deletion.append(('af', determineAlleleFrequency(phase))) 				
		deletion.append(('phase', phase))		
		deletion = dict(deletion)	
		if deletion['af'] <= majorAF_threshold and 1 - deletion['af'] <= majorAF_threshold:		
			deletions.append(deletion)

	vcf_file.close()
	return deletions 	

# returnDeletions returns only the data for deletions 
#def returnDeletions (vcf_filename, indices, majorAF_threshold = 1):
#	vcf_file = openVCFFile(vcf_filename)
#	discardVCFHeaders(vcf_file)
#	vcf_file.readline() # discard the column description header
#	
#	deletions = [] # data is stored here
#
	# loop over all lines until you find a deletion
#	for variant in vcf_file.readlines():
#		info = [value for value in variant.split()]
#		if len(info[3]) > 1: # this is a deletion	
#			deletion = [('pos', int(info[1])), ('length', len(info[3]) - 1)] 
#			#deletion.append(('af', float(info[7].split(';')[1].split('=')[1])))
#			phase = returnPhase(info, indices)
#			deletion.append(('maf', determineMajorAlleleFrequency(phase))) 				
#			deletion.append(('phase', phase))		
#			deletion = dict(deletion)	
#			if deletion['maf'] <= majorAF_threshold:		
#				deletions.append(deletion)				
#	vcf_file.close()
#	return deletions

def print2x2Table (x, y):
	n11, n12, n21, n22 = 0, 0, 0, 0 # Table 
	for i in range(len(x)):
		if x[i] == 1 and y[i] == 1:
			n11 = n11 + 1 
		elif x[i] == 1 and y[i] == 0:
			n21 = n21 + 1
		elif x[i] == 0 and y[i] == 1:
			n12 = n12 + 1
		else:
			n22 = n22 + 1  
	print '\tx=1\tx=0\t|\ttotal' 
	print '-----------------------------------------------'
	print 'y=1\t', n11, '\t', n12, '\t|\t', n11 + n12
	print 'y=0\t', n21, '\t', n22, '\t|\t', n21 + n22
	print '-----------------------------------------------'
	print '\t', n11 + n21, '\t', n12 + n22, '\t|\t', len(x)


def returnPositions (vcf_filename, indices, positions):
	vcf_file = openVCFFile(vcf_filename)
	discardVCFHeaders(vcf_file)
	vcf_file.readline() # discard the column description header

	variants = [] # data is stored here
	
	# loop over all lines and store all the variants on the given positions
	for variant in vcf_file.readlines():
		info = [value for value in variant.split()]
		position = int(info[1])
		if position in positions:
			phase = returnPhase(info, indices)
			af = determineAlleleFrequency(phase) 
			variants.append(dict([('pos', position), ('af', af), ('phase', phase)]))				

	return variants 
			

