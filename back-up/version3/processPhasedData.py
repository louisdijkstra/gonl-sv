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
# Date: December 2013
import math 
from scipy.stats.stats import pearsonr
import vcf as vcf 
import hg19 as hg19
import gamma as gamma

DIST_THRESHOLD 		= 2000000 
MAJOR_AF_THRESHOLD 	= .96

# Returns the distance between the SNP and the deletion (CNV)
# Returns 0 when the SNP lies within the deletion
def distanceSNP_CNV (snp_pos, cnv_start, cnv_length):
	if snp_pos < cnv_start:
		return cnv_start - snp_pos
	elif snp_pos > cnv_start + cnv_length:
		return snp_pos - cnv_start + cnv_length 
	else: # snp overlaps with the deletion
		return 0 

# Writes arrays to a prescripted output file 
def writeToOutputFile(output, chr):
	filename = 'chr' + str(chr) + '.txt' 
	output_file = open(filename, 'a')
	for line in output:
		for item in line:
			output_file.write(str(item) + '\t')
		output_file.write('\n')

	output_file.close() 


# returns the snp deletion pair with the highest R^2
# in addition returns the minimal distance from a SNP to a deletion
def returnSnpDeletionPair (snp_pos, snp_phase, deletions):
	index_del_max = -1
	max_r = 0.0 
	final_p = None
	dist_snp_del = None
	min_dist_snp_del = float('inf')

	for i in range(len(deletions)):
		deletion = deletions[i]
		current_dist_snp_del = distanceSNP_CNV(snp_pos, deletion['pos'], deletion['length'])
		if current_dist_snp_del <= DIST_THRESHOLD:
			if min_dist_snp_del > current_dist_snp_del:
				min_dist_snp_del = current_dist_snp_del 

			[R, P] = pearsonr(snp_phase, deletion['phase'])
			if not math.isnan(R) and abs(max_r) <= abs(R):
				max_r = R 
				final_p = P 
				dist_snp_del = current_dist_snp_del
				index_del_max = i 

	if index_del_max != -1:
		return [max_r, final_p, deletions[index_del_max], dist_snp_del, min_dist_snp_del] 		
	else:
		return [max_r, final_p, None, dist_snp_del, min_dist_snp_del]

# reads the hit SNPs from GWAS_catalogue_caucasian_SNPIDs.alleles.txt
# Hit SNPs for which the hit allele is not defined are discarded and stored in discarded_hit_snps_positions
def readHitSNPs (chromosome):
	snp_files = open('/ufs/dijkstra/Projects/SNPs_CNVs/data/GWAS_catalogue_caucasian_SNPIDs.alleles.txt', 'r')
	snp_files.readline() # discard header
	
	hit_snps_positions = [] # positions of the hit snps are stored here
	snps = [] # positions and hit allele are stored here
	discarded_hit_snps_positions = [] 

	chromosome = 'chr' + str(chromosome)

	for line in snp_files.readlines():
		info = [value for value in line.split()]
		if info[0] == chromosome: # snp lies on the chromosome of interest
			if len(info) == 4: # check whether the hit allele is given in the file
				hit_snps_positions.append(int(info[1]))
				snps.append( dict( [('pos', int(info[1])), ('hit_allele', info[3])] ) )
			else:
				discarded_hit_snps_positions.append(int(info[1]))
	
	snp_files.close()
	return [hit_snps_positions, snps, discarded_hit_snps_positions]


def processHitSNP (snp_pos, data, hit_snps, indices_individuals):
	for snp in hit_snps:
		if snp['pos'] == snp_pos:
			phase = vcf.returnPhaseWithAllele (data, indices_individuals, snp['hit_allele']) 
			snp_af = vcf.determineAlleleFrequency(phase)
			return [phase, snp_af]
		

# reads in the necessary data for preprocessing the data
def gatherData (vcf_filename, chromosome, individuals):
	[hit_snps_positions, hit_snps, discarded_snps_positions] = readHitSNPs (chromosome) 
	genes = hg19.read (chromosome)
	indices_individuals = vcf.returnColumns (vcf_filename, individuals) 
	deletions = vcf.returnDeletions (chromosome, indices_individuals, MAJOR_AF_THRESHOLD)
	return [hit_snps_positions, hit_snps, discarded_snps_positions, genes, indices_individuals, deletions]
	

	
def preprocessData (vcf_filename, chromosome, individuals=('parents', 'children')):
	[hit_snps_positions, hit_snps, discarded_snps_positions, genes, indices_individuals, deletions] = gatherData(vcf_filename, chromosome, individuals)	
	
	vcf_file = vcf.openVCFFile(vcf_filename)
	vcf.discardVCFHeaders(vcf_file)
	vcf_file.readline() # discard the column description header

	output = [] # results are stored here and later written to file
	n_hit_snps, n_discarded_snps, n_nonhit_snps = 0, 0, 0 

	for variant in vcf_file.readlines (): # move through the phased variants in this file
		data = [value for value in variant.split()]
		if len(data[3]) == 1: # variant is a SNP and no deletion
			snp_pos = int(data[1])
			snp_type = hg19.determineType(snp_pos, genes)
			dist_tss = hg19.distanceToTSS(snp_pos, genes)
			if snp_pos in hit_snps_positions: # SNP is a hit SNP
				[phase, snp_af] = processHitSNP (snp_pos, data, hit_snps, indices_individuals)
				[max_r, p, d, dist_snp_del, min_dist_snp_del] = returnSnpDeletionPair (snp_pos, phase, deletions)	
				if d is not None:
					[g,g_n] = gamma.returnGamma(phase, d['phase'], snp_af, d['af'])
					[a,b,c,e] = vcf.return2x2Table(phase, d['phase']) 
					n_hit_snps = n_hit_snps + 1
					output.append([chromosome, 1, snp_pos, snp_type, snp_af, d['pos'], d['length'], d['af'], max_r, p, dist_snp_del, dist_tss, min_dist_snp_del, g, g_n, a, b, c, e])		
				else:
					n_discarded_snps = n_discarded_snps + 1
					discarded_snps_positions.append(snp_pos) 
			elif snp_pos not in discarded_snps_positions: # non-hit SNP 
				phase 	= vcf.returnPhase(data, indices_individuals) 
				snp_af 	= vcf.determineAlleleFrequency (phase)
				[max_r, p, d, dist_snp_del, min_dist_snp_del] = returnSnpDeletionPair (snp_pos, phase, deletions)	
				if d is not None:
					[g,g_n] = gamma.returnGamma(phase, d['phase'], snp_af, d['af'])
					[a,b,c,e] = vcf.return2x2Table(phase, d['phase']) 
					n_nonhit_snps = n_nonhit_snps + 1
					output.append([chromosome, 0, snp_pos, snp_type, snp_af, d['pos'], d['length'], d['af'], max_r, p, dist_snp_del, dist_tss, min_dist_snp_del, g, g_n, a, b, c, e])		
				else:
					n_discarded_snps = n_discarded_snps + 1
					discarded_snps_positions.append(snp_pos)  
	writeToOutputFile(output, chromosome)
	vcf_file.close()	
	print '# hit SNPs: ', n_hit_snps, '\t#discarded SNPs:', n_discarded_snps, '\t# non hit SNPs: ', n_nonhit_snps

	
def main():
	print 'Chromosome 1'
	vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr1_snps_and_dels.vcf'
	preprocessData(vcf_filename, 1)
	for chr in range(2,22):
		print 'Chromosome ' + str(chr)
		vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + str(chr) + '_snps_and_dels.vcf.gz'
		preprocessData(vcf_filename, chr)
	print 'Chromosome 22'
	vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr22_snps_and_dels.vcf'
	preprocessData(vcf_filename, 22)

#main()
print 'Chromosome 22'
vcf_filename = '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr22_snps_and_dels.vcf'
preprocessData(vcf_filename, 22, individuals=('children'))

