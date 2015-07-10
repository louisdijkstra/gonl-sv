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

# processPhasedData.py processes the phased data of SNPs and deletions
#	from the GoNL project. Here, we focus only on the data of 
#	independent children (we consider only one child from each twin pair).
#	
#	For every SNP, we find the deletion which is on the same chromosome
#	and less than 2 million bp away , that results in the highest R^2. 
#	We store for every (SNP, deletion)-pair, their positions, length,
#	allele frequencies, Pearson's R and p-value. 
#	We restrict ourselves to autosomes (chromosome 1-22).
#
#	The results are stored in 
#
#		snp_del_maxR2_chr?.phased 	where ? is the number of the
#						chromosome
#
#	The columns in the output files represent:
#		#CHR 	- the chromosome on which the SNP and deletion pair lie
#		POS 	- 
# 
# Author: Louis Dijkstra
# Date: October 2013
import math 
from scipy.stats.stats import pearsonr
import vcf as vcf 
import hg19 as hg19

DIST_THRES = 2000000 
MAJOR_AF_THRESHOLD = .96

# Returns the distance between the SNP and the deletion (CNV)
# Returns -1 and prints a warning when the SNP lies within the deletion
def distanceSNP_CNV (snp_pos, cnv_start, cnv_length):
	if snp_pos < cnv_start:
		return cnv_start - snp_pos - 1
	elif snp_pos > cnv_start + cnv_length - 1:
		return snp_pos - cnv_start - cnv_length 
	else:
		print 'WARNING: SNP at position ', snp_pos, ' appears to have been deleted' 
		return -1 


# Writes one array to predescripted output file
def writeLineToOutputFile(output, chr):
	filename = 'NEW_snp_del_pearson_r_chr' + str(chr) + '.phased' 
	output_file = open(filename, 'a')
	for item in output:
		output_file.write(str(item) + '\t')
	output_file.write('\n')

	output_file.close() 

# Writes arrays to a prescripted output file 
def writeToOutputFile(output, chr):
	filename = 'androniki_r_chr' + str(chr) + '.phased' 
	output_file = open(filename, 'a')
	for line in output:
		for item in line:
			output_file.write(str(item) + '\t')
		output_file.write('\n')

	output_file.close() 

def maximumR_hitsnp (snp_pos, snp_phase, deletions):
	index_del_max = -1 
	max_r = 0.0 
	final_p = None

	for i in range(len(deletions)):
		deletion = deletions[i]
		dist_snp_del = distanceSNP_CNV(snp_pos, deletion['pos'], deletion['length'])
		if dist_snp_del <= DIST_THRES:
			[R, P] = pearsonr(snp_phase, deletion['phase']) 
			if not math.isnan(R) and max_r < R:
				max_r = R
				final_p = P 
				index_del_max = i 
	
	if index_del_max != -1: 			 
		return [max_r, final_p, deletions[index_del_max]]
	else:
		return [-2,-2,None]

def maximumR (snp_pos, snp_phase, deletions):
	index_del_max = -1 
	max_r = 0.0 
	final_p = None

	for i in range(len(deletions)):
		deletion = deletions[i]
		dist_snp_del = distanceSNP_CNV(snp_pos, deletion['pos'], deletion['length'])
		if dist_snp_del <= DIST_THRES:
			[R, P] = pearsonr(snp_phase, deletion['phase']) 
			if not math.isnan(R) and abs(max_r) < abs(R):
				max_r = R
				final_p = P 
				index_del_max = i 
	
	if index_del_max != -1: 			 
		return [max_r, final_p, deletions[index_del_max]]
	else:
		return [-2,-2,None]


def readHitSNPs (chromosome):
	snp_files = open('/ufs/dijkstra/Projects/SNPs_LD_deletions_phased/data/GWAS_catalogue_caucasian_SNPIDs.alleles.txt', 'r')
	snp_files.readline() # discard header
	
	hit_snp_positions = [] 
	snps = [] # data is stored here

	chromosome = 'chr' + str(chromosome)

	for line in snp_files.readlines():
		info = [value for value in line.split()]
		if info[0] == chromosome and len(info) == 4: # check whether this is the chromosome we are interested in
			snps.append( dict( [('pos', int(info[1])), ('hit_allele', info[3])] ) )
			hit_snp_positions.append(int(info[1]))
	
	snp_files.close()
	return [hit_snp_positions, snps]


def processHitSNP (snp_pos, data, hit_snps, indices_children):
	for snp in hit_snps:
		if snp['pos'] == snp_pos:
			phase = vcf.returnPhaseWithAllele (data, indices_children, snp['hit_allele'])
			snp_af = vcf.determineAlleleFrequency(phase)
			return [phase, snp_af]
		
	
	



def processPhasedData (vcf_filename, chromosome):
	[hit_snp_positions, hit_snps] = readHitSNPs(chromosome)
	genes = hg19.read(chromosome)
	indices_children = vcf.returnColumnsOfChildren(vcf_filename)	
	deletions = vcf.returnDeletions(chromosome, indices_children, MAJOR_AF_THRESHOLD)
	vcf_file = vcf.openVCFFile(vcf_filename)
	vcf.discardVCFHeaders(vcf_file)
	vcf_file.readline() # discard the column description header

	output = [] # output is stored here
	
	for variant in vcf_file.readlines():
		data = [value for value in variant.split()]
		if len(data[3]) == 1: # variant is a SNP
			# then determine whether this is a hit snp:
			snp_pos = int(data[1])
			#if snp_pos in hit_snp_positions: # hit snp
			#	[phase, snp_af] = processHitSNP (snp_pos, data, hit_snps, indices_children)
			#	if snp_af <= MAJOR_AF_THRESHOLD and 1 - snp_af <= MAJOR_AF_THRESHOLD:
			#		snp_type 	= hg19.determineType(snp_pos, genes)
			#		dist_tss 	= hg19.distanceToTSS(snp_pos, genes)		
			#		[max_r, p, d] 	= maximumR_hitsnp (snp_pos, phase, deletions)
			#		#print '1', max_r, p
			#		if d is not None:
			#			dist_snp_del = distanceSNP_CNV(snp_pos, d['pos'], d['length']) 
			#			output.append([chromosome, 1, snp_pos, snp_type, snp_af, d['pos'], d['length'], d['af'], max_r, p, dist_snp_del, dist_tss])
			#else: # non-hit snp
			phase 	= vcf.returnPhase(data, indices_children) 
			snp_af 	= vcf.determineAlleleFrequency (phase)
			if snp_af <= MAJOR_AF_THRESHOLD and 1 - snp_af <= MAJOR_AF_THRESHOLD:
				snp_pos 	= int(data[1])
				snp_type 	= hg19.determineType(snp_pos, genes)
				dist_tss 	= hg19.distanceToTSS(snp_pos, genes)		
				[max_r, p, d] 	= maximumR (snp_pos, phase, deletions)
				type_snp = 0
				if snp_pos in hit_snp_positions:
					type_snp = 1 
				#print '0', max_r, p
				if d is not None:
					dist_snp_del = distanceSNP_CNV(snp_pos, d['pos'], d['length']) 
					output.append([chromosome, type_snp, snp_pos, snp_type, snp_af, d['pos'], d['length'], d['af'], max_r, p, dist_snp_del, dist_tss])
	writeToOutputFile(output, chromosome)
	vcf_file.close()	

	
def main():
	#print 'Chromosome 1'
	#vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr1_snps_and_dels.vcf'
	#processPhasedData(vcf_filename, 1)
	for chr in range(9,22):
		print 'Chromosome ' + str(chr)
		vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chr) + '_snps_and_dels.vcf.gz'
		processPhasedData(vcf_filename, chr)


#main()
vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr20_snps_and_dels.vcf.gz'
processPhasedData(vcf_filename, 20)

