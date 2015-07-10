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
from scipy.stats.stats import pearsonr
import vcf as vcf 
import hg19 as hg19
import gamma as gamma
import processPhasedData as pr

DIST_THRESHOLD 		= 1000000
MAJOR_AF_THRESHOLD 	= .96 

# Returns the distance between the SNP and the deletion (CNV)
# Returns 0 when the SNP lies within the deletion
def distanceSNP_Deletion (snp_pos, d_start, d_length):
	if snp_pos < d_start:
		return d_start - snp_pos
	elif snp_pos > d_start + d_length:
		return snp_pos - d_start + d_length 
	else: # snp overlaps with the deletion
		return 0 

# Writes SNP data to a prescripted output file 
def writeSNPToOutputFile(snp_data, chromosome):
	print snp_data
	filename = 'all_SNP_Del_combinations_' + str(DIST_THRESHOLD) + 'chr' + str(chromosome) + '.txt' 
	output_file = open(filename, 'a')
	for item in snp_data[0]:
		output_file.write(str(item) + '\t')
	output_file.write('\n')
	for i in range(1,len(snp_data)):
		for item in snp_data[1]:
			output_file.write(str(item) + '\t')
		output_file.write('\n')

	output_file.close() 

def returnDeletions(snp_pos, snp_af, phase, deletions):
	del_output = [] 
	for d in deletions:
		d_snp_del = distanceSNP_Deletion (snp_pos, d['pos'], d['length'])
		if d_snp_del <= DIST_THRESHOLD and (1 - MAJOR_AF_THRESHOLD <= d['af'] <= MAJOR_AF_THRESHOLD):
			[R, P] = pearsonr(phase, d['phase'])
			[g, g_n] = gamma.returnGamma (phase, d['phase'], snp_af, d['af'])
			[a,b,c,e] = vcf.return2x2Table(phase, d['phase']) 
			del_output.append(['DEL', d['pos'], d['length'], d['af'], d_snp_del, R, g, g_n, a, b, c, e])
	return del_output

def processSNPs (vcf_filename, chromosome, individuals=('parents', 'children')):
	# gather the required data
	[hit_snps_positions, hit_snps, discarded_snps_positions] = pr.readHitSNPs (chromosome)
	genes = hg19.read (chromosome)
	indices_individuals = vcf.returnColumns (vcf_filename, individuals) 
	deletions = vcf.returnDeletions (chromosome, indices_individuals, MAJOR_AF_THRESHOLD)
	# open file	
	vcf_file = vcf.openVCFFile(vcf_filename)
	vcf.discardVCFHeaders(vcf_file)
	vcf_file.readline() # discard the column description header
	
	for variant in vcf_file.readlines(): # move through the genetic variants
		snp_output = [] 
		data = [value for value in variant.split()]
		if len(data[3]) == 1: # variant is a SNP and no deletion
			snp_pos = int(data[1])
			snp_type = hg19.determineType(snp_pos, genes)
			dist_tss = hg19.distanceToTSS(snp_pos, genes)	
			if snp_pos in hit_snps_positions: # SNP is a hit SNP
				[phase, snp_af] = pr.processHitSNP (snp_pos, data, hit_snps, indices_individuals)
				if 1 - MAJOR_AF_THRESHOLD <= snp_af <= MAJOR_AF_THRESHOLD:
					snp_output.append(['HITSNP', snp_pos, snp_af, snp_type, dist_tss]) 
					dels = returnDeletions(snp_pos, snp_af, phase, deletions) # TODO implement
					if len(dels) > 0:
						snp_output.append(dels)
			elif snp_pos not in discarded_snps_positions:
				phase 	= vcf.returnPhase(data, indices_individuals) 
				snp_af 	= vcf.determineAlleleFrequency (phase)
				if 1 - MAJOR_AF_THRESHOLD <= snp_af <= MAJOR_AF_THRESHOLD:				
					snp_output.append(['NONHITSNP', snp_pos, snp_af, snp_type, dist_tss]) 
					dels = returnDeletions(snp_pos, snp_af, phase, deletions) # TODO implement
					if len(dels) > 0:
						snp_output.append(dels)
			if len(snp_output) > 0:
				writeSNPToOutputFile(snp_output, chromosome)
	vcf_file.close()

print 'Chromosome 22'
vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr22_snps_and_dels.vcf'
processSNPs(vcf_filename, 22, individuals=('children'))

