/*
 * Copyright (C) 2015 Louis Dijkstra
 * 
 * This file is part of gonl-sv
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * tagsnps_main.c
 * 
 * Author: Louis Dijkstra
 * E-mail: dijkstra@cwi.nl
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tagsnps_input.h"
#include "haplotypes.h"

#define N_SNPS	1000 


size_t d(size_t pos1, size_t pos2) {
	/* return distance between two positions */
	return abs(pos1 - pos2) ; 
}

int main(int argc, char *argv[])
{
	// parse the command line   
	char input_filename1[256]; 
	char input_filename2[256]; 
	parameters p ; 
	parse_arguments(&p, input_filename1, input_filename2, argc, argv) ;
	
	if (p.verbose == 1) {
		print_parameters(&p) ;
	} 

	FILE *fp_del = fopen(input_filename1, "r") ; 
	FILE *fp_snp = fopen(input_filename2, "r") ; 
	
	if (fp_del == NULL) {
		printf("Could not open file %s\n", input_filename1) ; 
		exit(EXIT_FAILURE) ; 
	}
	if (fp_snp == NULL) {
		printf("Could not open file %s\n", input_filename2) ; 
		exit(EXIT_FAILURE) ; 
	}

	variant del ;  
	size_t length_snps = N_SNPS ; 
	variant * snps = malloc(length_snps * sizeof(variant)) ; 
	size_t i, discard, n_snps = 0 ; 
	size_t status_del, status_snp ; 

	int index_max ; 
	double max_r2 = 0.0 ; 
	size_t A,B,C,D;
	char allele1, allele2 ; 
	double p_del_allele1, p_del_allele2; 
	while(1) {
		// read in the data 
		status_del = obtainVariant(fp_del, &del) ; 
		if (status_del == END_OF_FILE_REACHED) {
			break ; 
		}
			
		discard = 0 ; 
		for (i = 0; i < n_snps; i ++) 
			if (d(del.position, snps[i].position) > p.search_range) 
				discard ++ ; 
			else 
				break ;
		
		if (discard != 0) {
			for (i = discard; i < n_snps; i ++) {
				snps[i - discard] = snps[i] ; 
			}
			n_snps -= discard ; 
		}

		while (n_snps == 0 || d(snps[n_snps - 1].position, del.position) < p.search_range) {
			variant snp ; 			
			status_snp  = obtainVariant(fp_snp, &snp) ;
			if (status_snp == END_OF_FILE_REACHED) 
				break ; 
			snps[n_snps] = snp ; 
			n_snps ++ ; 
			if (n_snps == length_snps) {
				length_snps += N_SNPS ; 
				snps = realloc(snps, sizeof(variant) * length_snps);
				if (snps == NULL) {
					printf("ERROR: insufficient memory for reallocation.") ; 
					exit(EXIT_FAILURE) ; 
				}		
			}
		}
		
		double * p_values = malloc(n_snps * sizeof(double)) ; 
		double * R = malloc(n_snps * sizeof(double)) ; 
		if (p_values == NULL || R == NULL) {
			printf("ERROR: insufficient memory.\n") ; 
			exit(EXIT_FAILURE) ; 
		}
	
		for (i = 0; i < n_snps; i ++) {
			R[i] = haplo_pearsonsR(del.haplotypes, del.af, snps[i].haplotypes, snps[i].af, del.n) ;
			p_values[i] = haplo_fishersTest(del.haplotypes, snps[i].haplotypes, del.n) ; 
		}	
		
		index_max = -1 ; 
		max_r2 = 0 ; 

		for (i = 0; i < n_snps; i ++) {
			if (pow(R[i], 2) >= p.r2_threshold && p_values[i] <= p.p_threshold) { 
				if (pow(R[i], 2) > max_r2) {
					index_max = i ;
					max_r2 = pow(R[i], 2) ; 
				}
			}
		}



		if (index_max != -1) {
			printf("%zd\t%zd\t%zd\t", del.autosome, del.position, del.length) ;
			if (snps[index_max].af < .5) {
				snps[index_max].af = 1.0 - snps[index_max].af ; 
				for (i = 0; i < snps[index_max].n; i ++) {
					if (snps[index_max].haplotypes[i] == 1) 
						snps[index_max].haplotypes[i] = 0 ;
					else 
						snps[index_max].haplotypes[i] = 1 ;
				}
				allele1 = snps[index_max].alt ; 
				allele2 = snps[index_max].ref ; 
				R[index_max] *= -1.0 ; 
			} else {
				allele1 = snps[index_max].ref ; 
				allele2 = snps[index_max].alt ; 
			}
			printf("%zd\t%c\t%c\t%lf\t%e\t", snps[index_max].position, allele1, allele2, R[index_max], p_values[index_max]) ; 
			haplo_2x2table (snps[index_max].haplotypes, del.haplotypes, del.n, &A, &B, &C, &D) ;
			
			p_del_allele1 = (double)A / (A + C) ; 
			p_del_allele2 = (double)B / (B + D) ; 				
			printf("%lf\t%lf\t%zd\t%zd\t%zd\t%zd\n", p_del_allele1, p_del_allele2, A, B, C, D) ; 
		} 
	}
	
    	fclose(fp_del);
	fclose(fp_snp);
    	exit(EXIT_SUCCESS) ; 
}
