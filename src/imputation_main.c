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
 * imputation_main.c
 * 
 * Author: Louis Dijkstra
 * E-mail: dijkstra@cwi.nl
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "imputation_input.h"
#include "haplotypes.h"

#define N_SNPS	1000

size_t d(size_t pos1, size_t pos2) {
	/* return distance between two positions */
	return abs(pos1 - pos2) ; 
}

double logit(double x) {
	return log(x) - log(1.0 - x) ; 
}

double perform_leave_one_out (variant* deletion, variant* snps, size_t n, parameters* p) {
	if (n > 0 && d(snps[n-1].position, deletion->position) > p->search_range) 
		n -- ; 

	if (n == 0)  /* no data to go on */
		return 0 ; 
	
	size_t i, j, k ; 
	size_t m = deletion->n ; // number of observations

	double m0 = 0 ; 
	for (j = 0; j < m; j++) 
		if (deletion->haplotypes[j] == 0) 
			m0 ++ ; 
	double m1 = m - m0 ; 
	double eta = m1 / m, eta_k ; 
	size_t a=0, b=0, c=0, d=0 ; 
	size_t number_correct = 0 ; 
	if (eta == 0.0) {
		// print all calls to be zero
		for (k = 0; k < m; k ++) {
			if (deletion->haplotypes[k] == 1) {
				c ++ ; 
			} else {
				d ++ ; 
				number_correct ++ ; 
			}
		}
		printf("- %zd %zd %zd ", deletion->autosome, deletion->position, deletion->length) ; 
		printf("%zd %zd %zd %zd %zd ", n, a, b, c, d) ; 
		printf("%f\n", (double) number_correct / m * 100) ;
		return (double) number_correct / m * 100 ; 
	} else if (eta == 1.0) {
		// print all calls to be one
		for (k = 0; k < m; k ++) {
			if (deletion->haplotypes[k] == 1) {
				a ++ ; 
				number_correct ++ ; 
			} else {
				b ++ ; 
			}
		}
		printf("- %zd %zd %zd ", deletion->autosome, deletion->position, deletion->length) ; 
		printf("%zd %zd %zd %zd %zd ", n, a, b, c, d) ; 
		printf("%f\n", (double) number_correct / m * 100) ;
		return (double) number_correct / m * 100 ; 
	}

	/* allocate memory for the parameters */
	double * pi0 = calloc(n, sizeof(double)) ; 
	double * pi1 = calloc(n, sizeof(double)) ; 
	double * pi0_k = malloc(n*sizeof(double)) ; 
	double * pi1_k = malloc(n*sizeof(double)) ; 
	if (pi0 == NULL || pi1 == NULL || pi0_k == NULL || pi1_k == NULL) {
		printf("ERROR: insufficient memory.") ; 
		exit(EXIT_FAILURE) ; 
	}

	for (i = 0; i < n; i ++) {
		for (j = 0; j < m; j ++) {
			if (deletion->haplotypes[j] == 1) 
				pi1[i] += snps[i].haplotypes[j] ; 
			else 
				pi0[i] += snps[i].haplotypes[j] ; 
		}
	}

	for (i = 0; i < n; i ++) {
		pi0[i] /= m0 ; 
		pi1[i] /= m1 ; 
	}

	/* allocate memory for classification */
	double sufficient_statistic = 0.0 ;  
	double threshold = 0.0 ; 
	size_t call ; 
	for (k = 0; k < m; k ++) {
		eta_k = eta + (1.0 / (m - 1)) * (eta - deletion->haplotypes[k]) ; 
		if (eta_k == 0.0) {
			call = 0 ; 
		} else if (eta_k == 1.0) {
			call = 1 ;
		} else {
			
			for (i = 0; i < n; i ++) {
				if (deletion->haplotypes[k] == 1) {
					pi0_k[i] = pi0[i] ; 
					pi1_k[i] = (1 / (m1 - 1)) * (m1 * pi1[i] - snps[i].haplotypes[k]) ; 
				} else {
					pi0_k[i] = (1 / (m0 - 1)) * (m0 * pi0[i] - snps[i].haplotypes[k]) ; 
					pi1_k[i] = pi1[i] ;  
				}
			}			
			
			sufficient_statistic = 0.0 ;
			threshold = -1 * logit(eta_k) ; 
			call = 2 ; 
			for (i = 0; i < n; i ++) {
				if (snps[i].haplotypes[k] == 1) {
					if (pi1_k[i] == 1) {
						// must be 1
						call = 1 ; 
						break ; 
					}
					if (pi1_k[i] == 0) {
						call = 0 ; 
						// must be 0
						break ; 
					}
					if (pi0_k[i] == 1) {
						call = 0 ; 
						// must be 0
						break ; 
					}
					if (pi0_k[i] == 0) {
						call = 1 ; 
						// must be 1
						break ; 
					}
					sufficient_statistic += (logit(pi1_k[i]) - logit(pi0_k[i])) ; 	
				}
				if (pi1_k[i] == 1.0) {
					call = 1 ; 
					break ; 
				}
				if (pi0_k[i] == 1.0) {
					call = 0 ; 
					break ; 
				}
				threshold += log(1.0 - pi0_k[i]) - log(1.0 - pi1_k[i]) ; 
			}	
			
			if (call == 2) {
				if (sufficient_statistic >= threshold) {
					call = 1 ; 
				} else {
					call = 0 ; 
				}
			}

		}

		if (deletion->haplotypes[k] == 1) {
			if (call == 1) 
				a ++ ; 
			else
				b ++ ; 
		} else {
			if (call == 1) 
				c ++ ; 
			else
				d ++ ; 
		}

		if (call == deletion->haplotypes[k]) 
			number_correct ++ ; 
	}

	printf("- %zd %zd %zd ", deletion->autosome, deletion->position, deletion->length) ; 
	printf("%zd %zd %zd %zd %zd ", n, a, b, c, d) ; 
	printf("%f\n", (double) number_correct / m * 100) ; 

	return (double) number_correct / m * 100 ; 

} 


int main(int argc, char *argv[])
{
	/* parse the command line */ 
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

	variant * subset_snps = malloc(p.number_of_snps * 2 * sizeof(variant)) ; 
	size_t n_subset_snps = 0 ; 

	size_t i, discard, n_snps = 0 ; 
	size_t status_del, status_snp ; 

	size_t n_deletions = 0 ; 
	double percentage_correct = 0.0 ; 

	while(1) {
		/* read in the deletion */
		status_del = obtainVariant(fp_del, &del) ; 
		if (status_del == END_OF_FILE_REACHED) {
			break ; 
		}

		/* read in the SNPs */
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
		
		n_deletions ++ ; 
		if (p.number_of_snps == 0) 
			percentage_correct += perform_leave_one_out (&del, snps, n_snps, &p) ; 
		else {
			int index_start = 0, index_end = 0 ; 
			for (i = 0; i < n_snps; i ++) {
				//printf("snp position %zd\t deletion %zd\n", snps[i].position, del.position) ; 
				if (snps[i].position > del.position) {
					index_start = i - p.number_of_snps ; 
					if (index_start < 0)
						index_start = 0 ; 
					index_end = i + p.number_of_snps ; 
					if (index_end > n_snps) 
						index_end = n_snps ;  
					break ; 
				}
			}

			n_subset_snps = 0 ; 
			for (i = index_start; i < index_end ; i ++ ) {
				variant copy_snp = snps[i] ; 
				subset_snps[i - index_start] = copy_snp ; 
				n_subset_snps ++ ; 
			}
			//printf("%zd %zd %zd\n", index_start, index_end, n_subset_snps) ; 
			percentage_correct += perform_leave_one_out (&del, subset_snps, n_subset_snps, &p) ;
		}
		

		//printf("> total number of deletions: %zd\tpercentage correct: %f\n", n_deletions, percentage_correct / n_deletions) ; 

		/*
		printf("- %zd %zd %zd\n", del.autosome, del.position, del.length) ; 

		if (n_snps == 0) { // there is no data to go on
			printf("0\n") ; 
			printf(".\n") ; 
			for (i = 0; i < del.n; i ++) 
				printf("%zd:.:. ", del.haplotypes[i]) ; 	
		} else {
			//printf("%d ", n_snps) ;	
			//for (i = 0; i < n_snps-1; i ++) 
			//	printf("%zd ", snps[i].position) ; 
			//printf("%zd\n", snps[n_snps-1].position) ;

			
			double * pi0 = malloc(n_snps * sizeof(double)) ; 
			double * pi1 = malloc(n_snps * sizeof(double)) ; 
			int * m0 = malloc(n_snps * sizeof(int)) ; 
			int * m1 = malloc(n_snps * sizeof(int)) ; 

			double * pi0_a = malloc(n_snps * sizeof(double)) ; 
			double * pi1_a = malloc(n_snps * sizeof(double)) ; 

			if (pi0 == NULL || pi1 == NULL || m0 == NULL || m1 == NULL || pi0_a == NULL || pi1_a == NULL) {
				printf("ERROR: insufficient memory.") ; 
				exit(EXIT_FAILURE) ; 
			}
			
			
			for (i = 0; i < n_snps ; i ++) {
				m0[i] = 0 ; 
				m1[i] = 0 ; 
				for (j = 0; j < del.n; j ++) {
					if (del.haplotypes[j] == 1) {
						m1[i] ++ ; 
						pi1[i] += snps[i].haplotypes[j] ; 
					} else {
						m0[i] ++ ; 
						pi0[i] += snps[i].haplotypes[j] ; 
					}
				}
				
				if (m1[i] == 0) {
					pi1[i] = 0.0 ; 
				} else {
					pi1[i] /= m1[i] ; 
				} 

				if (m0[i] == 0) {
					pi0[i] = 0.0 ; 
				} else {
					pi0[i] /= m0[i] ; 
				} 
			}
			
		
			size_t k ; // 
			size_t m = del.n ; 
			double n_correct =0.0 ; 			
			for (k = 0; k < m; k++) {
				//af_est = del.af + (1.0 / (m - 1.0)) * (del.af - del.haplotypes[k]) ; // adjust af estimate
				if (del.haplotypes[k] == 0) {
					af_est *= m/(m-1) ; // adjust allele freq. estimate
					for (i = 0; i < n_snps; i ++) {
						pi1_a[i] = pi1[i] ; 
						pi0_a[i] = pi0[i] + (1.0 / (m0[i] - 1)) * (pi0[i] - snps[i].haplotypes[k]) ;
					}
				} else {
					af_est = (af_est*m - 1)/ (m-1) ; // adjust allele freq. estimate
					for (i = 0; i < n_snps; i ++) {
						pi0_a[i] = pi0[i] ;
						pi1_a[i] = pi1[i] + (1.0 / (m1[i] - 1)) * (pi1[i] - snps[i].haplotypes[k]) ;
					}
				}
				
				
				if (af_est == 0.0) {
					//printf("%zd:0:1.0 ", del.haplotypes[k]) ; 
				}
				double sum = 0.0 ; 
				int x ; 
				for (i = 0; i < n_snps; i ++) {
					x = snps[i].haplotypes[k] ; 
					if (x == 0 && pi1[i] == 1) {
						//printf("%zd:0:1.0 ", del.haplotypes[k]) ; 
						break ; 
					}
					if (x == 1 && pi1[i] == 0) {
						//printf("%zd:0:1.0 ", del.haplotypes[k]) ; 
						break ; 
					}
					sum += x * log(pi0_a[i] / pi1_a[i]) + (1 - x) * log((1-pi0_a[i]) / (1 - pi1_a[i])) ; 

				} 
				double post_prob = 1 / (1 + (1 - af_est)/af_est * exp(sum)) ; 
				if (post_prob >= 0.5) 
					//printf("%zd:1:%f ", del.haplotypes[k], post_prob) ; 
					if (del.haplotypes[k] == 1)
						n_correct ++ ; 
				else 
					//printf("%zd:0:%f ", del.haplotypes[k], 1.0 - post_prob) ; 
					if (del.haplotypes[k] == 0)
						n_correct ++ ;
				
			}
			printf("perc. correct %f\n", n_correct / del.n * 100.0) ; 
			
		}
		*/

	
		 
		

	}
	printf("> total number of deletions: %zd\tpercentage correct: %f\n", n_deletions, percentage_correct / n_deletions) ; 

	fclose(fp_del);
	fclose(fp_snp);
    	exit(EXIT_SUCCESS) ; 
}

		
		
