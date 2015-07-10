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
 * haplotypes.c
 * 
 * Author: Louis Dijkstra
 * E-mail: dijkstra@cwi.nl
 */
#include "haplotypes.h"

void haplo_print (size_t * haplotypes, size_t n) {
	/* Prints haplotypes to screen as a list */
	size_t i ;	
	printf("[") ; 
	for (i = 0 ; i < n - 1; i ++) {
		printf("%zd, ", haplotypes[i]) ; 
	}
	if (n != 0) {
		printf("%zd]\n", haplotypes[n-1]) ; 
	}
}	

double haplo_returnFrequency (size_t * haplotypes, size_t n) {
	/* Returns the frequency of ones in the binary list */
	size_t i ; 
	double ones = 0.0 ; 
	for (i = 0; i < n; i ++) 
		ones += haplotypes[i] ; 
	return ones / n ; 
}

double haplo_pearsonsR (size_t * haplotypes1, double mean1, size_t * haplotypes2, double mean2, size_t n) {
	/* Computes Pearson's R - returns 99 when R is not defined */
	size_t i ; 
	double std1 = 0.0, std2 = 0.0, r = 0.0 ;
	
	for (i = 0; i < n; i ++) {
		r += (haplotypes1[i] - mean1) * (haplotypes2[i] - mean2) ; 
		std1 += pow((haplotypes1[i] - mean1),2) ; 
		std2 += pow((haplotypes2[i] - mean2),2) ; 
	}
	std1 = sqrt(std1 / (n-1)) ; 
	std2 = sqrt(std2 / (n-1)) ; 
	if (std1 == 0.0 || std2 == 0.0) {
		return 99.0 ; 
	} 
	return r / (std1 * std2 * (n-1)) ; 
}

void haplo_2x2table (size_t * haplotypes1, size_t * haplotypes2, size_t n, size_t* a, size_t* b, size_t* c, size_t* d) {
	/* Constructs the 2x2 contigency table structured as follows: 
	 * 
 	 * 				haplotype2
	 * 			|	1	0	|  total
	 * 		-----------------------------------------
	 * 		1	|	a	c	|   a+c
	 * haplotype1	0	|	b	d	|   b+d 
	 * 		-----------------------------------------
	 *		total	| 	a+b	c+d	| a+b+c+d
	 */
	size_t i ; 
	
	// initialize
	(*a) = 0 ; 
	(*b) = 0 ; 
	(*c) = 0 ; 
	(*d) = 0 ; 

	for (i = 0; i < n; i ++) {
		if (haplotypes1[i] == 1) {
			if (haplotypes2[i] == 1)	(*a) += 1 ; 
			else				(*c) += 1 ; 
		} else {
			if (haplotypes2[i] == 1)	(*b) += 1 ; 
			else				(*d) += 1 ; 
		}
	}
}


double haplo_fishersTest (size_t * haplotypes1, size_t * haplotypes2, size_t n) {
	// construct table
	size_t a, b, c, d, first_row, first_column, k; 
	haplo_2x2table(haplotypes1, haplotypes2, n, &a, &b, &c, &d) ; 
	first_row = a + c ; 
	first_column = a + b ; 

	// probability of observing the actual table: 
	double p_obs_table = gsl_ran_hypergeometric_pdf(a, first_row, n - first_row, first_column) ; 

	double p_value = p_obs_table; 
	double p_table; 
	
	// walk through all possible tables 
	int min_a 	= first_row + first_column - n ; 
	size_t max_a 	= first_row ; 
	if (min_a < 0) 
		min_a = 0 ; 
	if (max_a > first_column) 
		max_a = first_column ; 
	
	for (k = min_a ; k < max_a ; k ++) {
		if (k != a) {
			p_table = gsl_ran_hypergeometric_pdf(k, first_row, n - first_row, first_column) ; 
			if (p_table <= p_obs_table) {
				p_value += p_table ; 
			}
		}
	}

	return p_value ; 
}

