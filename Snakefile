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

import os 
import os.path
import sys
import subprocess
import snakemake.utils

___author___ = "Louis Dijkstra"


gonl_deletions = "/data1/gonl/new-deletion-phasing/deletions-20-10000-phased-20140925.vcf"
human_omni_strand_file = "data/HumanOmni25-8v1-1_A-b37.strand"	
gwas_catalogue = "data/GWAS_catalogue_caucasian_SNPIDs.alleles.txt"
hg19_file = "data/hg19.txt"

### Parameters ###
MAF_THRESHOLD = ['0.04']
AUTOSOME = [str(i) for i in range(1,23)]
R2_THRESHOLD = ['0.8']  
P_THRESHOLD = ['0.05']

""" Program location """
SAMTOOLS 	= "/export/data1/gonl/bin/samtools"
BGZIP 		= "~/Projects/Software/tabix-0.2.6/bgzip"
TABIX 		= "~/Projects/Software/tabix-0.2.6/tabix" 


def GoNLSNPFileLocation(wildcards):
	if wildcards.autosome != '1' and wildcards.autosome != '22':
		return '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + wildcards.autosome + '_snps_and_dels.vcf.gz'
	return '/data2/as/hitSNP-CNV/data/phased-snps-dels/chr' + wildcards.autosome + '_snps_and_dels.vcf'
"""
rule annotate_gwas_snps:
	input:
		[ "data/gonl-snps/{name}.haplotypes" 
		, "data/gwas_catalogue.txt" ] 
	output:
		"data/gonl-snps/{name}.annotated.haplotypes"
	message:
		"Annotating the GWAS SNPs from {input[0]} with the SNPs in {input[1]}. Output is stored @ {output}."
	shell:
		"python bin/annotate-gwas-snps.py {input[0]} {input[1]} > {output}"
"""

rule create_all_pairs:
	input: 
		expand("/export/scratch3/dijkstra/gonl-sv/results/snp-deletion-pairs/snp_deletion.parents.chr{autosome}.a{maf_thres1}_{maf_thres2}.pairs",
			autosome = AUTOSOME,
			maf_thres1 = MAF_THRESHOLD,
			maf_thres2 = MAF_THRESHOLD)
	message:
		"Creating all SNP-deletion pairs for all autosomes."

rule create_pairs:
	input: 
		[ "data/gonl-snps/snps.parents.chr{autosome}.a{maf_thres1,[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?}.filtered.gwas.hg19.haplotypes"
		, "data/gonl-deletions/deletions.parents.chr{autosome}.a{maf_thres2,[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?}.haplotypes" ] 
	output:
		"/export/scratch3/dijkstra/gonl-sv/results/snp-deletion-pairs/snp_deletion.parents.chr{autosome}.a{maf_thres1}_{maf_thres2}.pairs"
	message:
		"Creating all SNP-deletion pairs given {input[0]} and {input[1]}. Output is stored @ {output}."
	shell:
		"src/gonl_create_pairs {input[0]} {input[1]} > {output}"
			
				

rule tag_all_deletions:
	input:
		expand("results/tag-snps/deletions.parents.chr{autosome}.a{maf_thres1}_{maf_thres2}.R2{r2_thres}.p{p_thres}.tagsnps", 
			autosome = AUTOSOME,
			maf_thres1 = MAF_THRESHOLD, 
			maf_thres2 = MAF_THRESHOLD,
			r2_thres = R2_THRESHOLD,
			p_thres = P_THRESHOLD) 
	message:
		"Tagging all GoNL Deletions."

rule tag_deletions:
	input:
		[ "data/gonl-deletions/deletions.parents.chr{autosome}.a{maf_thres1,[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?}.haplotypes"
		, "data/gonl-snps/snps.parents.chr{autosome}.a{maf_thres2,[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?}.filtered.gwas.hg19.haplotypes"]
	output:
		"results/tag-snps/deletions.parents.chr{autosome}.a{maf_thres1}_{maf_thres2}.R2{r2_thres}.p{p_thres}.tagsnps" 
	message:
		"Finding appropriate tag SNPs for the deletions on autosome {wildcards.autosome}."
	shell:
		"src/gonl_tag_deletions -r {wildcards.r2_thres} -p {wildcards.p_thres} {input[0]} {input[1]} > {output}"


rule annotate_hg19_all:
	input:
		expand("data/gonl-snps/snps.parents.chr{autosome}.a{maf_thres}.filtered.gwas.hg19.haplotypes", 
				autosome=AUTOSOME, 
				maf_thres=MAF_THRESHOLD)
	message:
		"Annotate all SNPs as exonic, intronic or intergenic."


rule annotate_hg19:
	input:
		[ "data/gonl-snps/{name}.gwas.haplotypes" 
		, hg19_file ] 
	output:
		"data/gonl-snps/{name}.gwas.hg19.haplotypes"
	message:
		"Determine type (intergenic/intronic/exonic) and the distance to the nearest TSS for every SNP in {input[0]} given {input[1]}. Output is stored @ {output}."
	shell:	
		"python bin/annotate-snps-type-tss.py {input[0]} {input[1]} > {output}"

rule split_deletions:
	input:
		expand("data/gonl-deletions/deletions.parents.chr{autosome}.a{maf_thres}.haplotypes",
				autosome=AUTOSOME,
				maf_thres=MAF_THRESHOLD) 
	message: "Split the deletion files into seperate files, one for each autosome."
	

rule filter_deletions_on_autosome:
	input:
		"data/gonl-deletions/deletions.parents.{specs}haplotypes"
	output:
		"data/gonl-deletions/deletions.parents.chr{autosome,\d+}.{specs}haplotypes"	
	message:
		"Filter out {input} all variants from autosome {wildcards.autosome}. Output is stored @ {output}."
	shell:
		"python bin/filter-haplotypes-on-autosome.py {input} {wildcards.autosome} > {output}"

rule annotate_all:
	input: 
		expand("data/gonl-snps/snps.parents.chr{autosome}.a{maf_thres}.filtered.gwas.haplotypes", 
				autosome=AUTOSOME, 
				maf_thres=MAF_THRESHOLD)
	message: 
		"Annotating all GoNL SNPs as being a GWAS or a non-GWAS SNP."

rule annotate_gwas_snps:
	input:
		[ "data/gonl-snps/{name}.haplotypes"
		, "data/gwas_catalogue.txt" ] 
	output:
		temp("data/gonl-snps/{name}.gwas.haplotypes")
	message:
		"Annotating the GWAS SNPs from {input[0]} with the SNPs in {input[1]}. Output is stored @ {output}."
	shell:
		"python bin/annotate-gwas-snps.py {input[0]} {input[1]} > {output}" 

rule filter_for_array_positions:
	input: 
		[ "data/{folder}/{name}.haplotypes" 
		, "data/snp_data.array" ] 
	output:
		temp("data/{folder}/{name}.filtered.haplotypes")
	message:
		"Filters out any position that is not on the array. Output is stored @ {output}."
	shell:
		"python bin/filter-haplotypes-on-position.py {input[0]} {input[1]} > {output}"
		

rule get_non_hit_snps:
	input: 
		[ "data/gonl-snps/{name}.haplotypes"
		, "data/gwas_catalogue.txt"]
	output:
		"data/gonl-snps/{name}.nonhit.haplotypes"
	message:
		"Extracting the non-hit SNPs from {input[0]}. Output is stored @ {output}."
	shell:
		"python bin/filter-haplotypes-on-position.py --leave-out {input[0]} {input[1]} > {output}"

rule filter_all_snp_haplotype_files:
	input:
		expand("data/gonl-snps/snps.parents.chr{autosome}.a{maf_thres}.filtered.haplotypes", 
				autosome=AUTOSOME, 
				maf_thres=MAF_THRESHOLD)
	message:
		"Filtering all the haplotype files for the SNPs."

rule create_snp_haplotype_files:	
	input: 
		expand("data/gonl-snps/snps.parents.chr{autosome}.a{maf_thres}.haplotypes", 
				autosome=AUTOSOME, 
				maf_thres=MAF_THRESHOLD)
	message:
		"Create all haplotype files for the SNPs."
	

rule extract_haplotypes_snps:
	output: temp("data/gonl-snps/snps.parents.chr{autosome,\d+}.a{maf_thres,[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?}.haplotypes")
	input: GoNLSNPFileLocation
	message: "Extracting the SNP haplotpyes of the parents for autosome {wildcards.autosome}. SNPs with a minor allele frequency below {wildcards.maf_thres} are discarded. Output is stored @ {output}."
	shell: "python bin/extract-haplotypes.py --parents --snps -a {wildcards.maf_thres} {input} > {output}"

rule extract_haplotypes_deletions:
	input: gonl_deletions
	output: "data/gonl-deletions/deletions.parents.haplotypes"
	message: "Extract the haplotypes of the parents."
	shell: "python bin/extract-haplotypes.py --parents --indels {input} > {output}"


### PREPROCESSING STEPS ###

rule preprocess_gwas_catalogue:
	input: gwas_catalogue
	output: "data/gwas_catalogue.txt"
	message: "Preprocessing and sorting the GWAS catalogue from {input}. Output is stored @ {output}."
	shell: "python bin/sort-gwas-catalogue.py {input} > {output}"

rule preprocess_strand_file:
	input: human_omni_strand_file 
	output: "data/snp_data.array"
	message: "Extracting and sorting the relevant information from {input}. Output is stored @ {output}."
	shell: "python bin/sort-strand-file.py {input} > {output}"



