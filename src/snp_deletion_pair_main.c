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

#include "snp_deletion_pair_input.h"
#include "haplotypes.h"

#define N_DELETIONS	100


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

	FILE *fp_snp = fopen(input_filename1, "r") ; 
	FILE *fp_del = fopen(input_filename2, "r") ; 
	
	if (fp_del == NULL) {
		printf("Could not open file %s\n", input_filename1) ; 
		exit(EXIT_FAILURE) ; 
	}
	if (fp_snp == NULL) {
		printf("Could not open file %s\n", input_filename2) ; 
		exit(EXIT_FAILURE) ; 
	}


	variant snp ; 
	size_t length_deletions = N_DELETIONS ; 
	variant * deletions = malloc(length_deletions * sizeof(variant)) ; 
	size_t i, discard, n_deletions = 0 ; 
	size_t status_snp, status_del ; 

	double R, p_value ; 
	size_t A, B, C, D ; 

	while(1) {
		// read in the data 
		status_snp = obtainVariant(fp_snp, &snp) ; 
		if (status_snp == END_OF_FILE_REACHED) {
			break ; 
		}

		// print data on SNP first
		printf("%zd %c %zd %c %c %c %zd ", 
				snp.autosome,
				snp.type, 
				snp.position,
				snp.ref,
				snp.alt,
				snp.hit_allele,
				snp.dist_tss) ;
		if (snp.functional == INTERGENIC) 
			printf("%s\n", INTERGENIC_STR) ; 
		else if (snp.functional == INTRONIC) 
			printf("%s\n", INTRONIC_STR) ;
		else 
			printf("%s\n", EXONIC_STR) ;
			
		discard = 0 ; 
		for (i = 0; i < n_deletions; i ++) 
			if (d(snp.position, deletions[i].position) > p.search_range) 
				discard ++ ; 
			else 
				break ;
		
		if (discard != 0) {
			for (i = discard; i < n_deletions; i ++) {
				deletions[i - discard] = deletions[i] ; 
			}
			n_deletions -= discard ; 
		}

		while (n_deletions == 0 || d(deletions[n_deletions - 1].position, snp.position) < p.search_range) {
			variant deletion ; 			
			status_del  = obtainVariant(fp_del, &deletion) ;
			if (status_del == END_OF_FILE_REACHED) 
				break ; 
			deletions[n_deletions] = deletion ; 
			n_deletions ++ ;

			// deletions 
			if (n_deletions == length_deletions) {
				length_deletions += N_DELETIONS ; 
				deletions = realloc(deletions, sizeof(variant) * length_deletions);
				if (deletions == NULL) {
					printf("ERROR: insufficient memory for reallocation.") ; 
					exit(EXIT_FAILURE) ; 
				}		
			}
		}


		
		for (i = 0; i < n_deletions; i ++) {
			if (d(deletions[i].position, snp.position) < p.search_range) {
				R = haplo_pearsonsR(deletions[i].haplotypes, deletions[i].af, snp.haplotypes, snp.af, snp.n) ;
				p_value = haplo_fishersTest(deletions[i].haplotypes, snp.haplotypes, snp.n) ;		
				haplo_2x2table (snp.haplotypes, deletions[i].haplotypes, snp.n, &A, &B, &C, &D) ;
				printf("- %zd %zd ", deletions[i].position, deletions[i].length) ; 
				printf("%lf %e %zd %zd %zd %zd\n", R, p_value, A, B, C, D) ; 
			}
		}	
	}
	
    	fclose(fp_snp);
	fclose(fp_del);
    	exit(EXIT_SUCCESS) ; 
}
