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

#ifndef HAPLOTYPES_H_
#define HAPLOTYPES_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_randist.h>

/* Prints haplotypes to screen */
void haplo_print (size_t * haplotypes, size_t n) ;

/* Returns the frequency of ones in the binary list */
double haplo_returnFrequency (size_t * haplotypes, size_t n) ; 

/* Computes Pearson's R - returns 99 when R is not defined */
double haplo_pearsonsR (size_t * haplotypes1, double mean1, size_t * haplotypes2, double mean2, size_t n) ; 

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
void haplo_2x2table (size_t * haplotypes1, size_t * haplotypes2, size_t n, size_t* A, size_t* B, size_t* C, size_t* D) ;

/* Performs Fishers Exact test */
double haplo_fishersTest (size_t * haplotypes1, size_t * haplotypes2, size_t n) ; 

#endif


