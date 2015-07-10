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

#ifndef SNP_DELETION_PAIR_INPUT_H_
#define SNP_DELETION_PAIR_INPUT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <unistd.h>
#include <float.h>

#define DELIM	" "	/* input file is space-delimited */

#define END_OF_FILE_REACHED 	1
#define VALID_VARIANT		0

// types of SNPs one can encouter:
#define INTERGENIC	0
#define INTRONIC	1
#define EXONIC		2

#define INTERGENIC_STR	"intergenic"
#define INTRONIC_STR	"intronic"	
#define EXONIC_STR	"exonic"

typedef struct {
	size_t search_range ; 
	short verbose ; 
} parameters ; 


/* 
 * The data structure (matches with one line)
 */ 
typedef struct {
	char type ; 		// '+' in case of an insertion, '-' for a deletion, '*' for a SNP and '!' for a GWAS SNP.
	size_t autosome ; 	// autosome that harbours the variant
	size_t position ; 	// position of the variant 	
	size_t length ; 	// length (in case of an indel)
	size_t n ; 		// number of haplotypes
	size_t * haplotypes ; 	// n binary observations
	double af ; 		// allele frequency
	char ref ; 
	char alt ; 

	// only used for SNPs
	char hit_allele ; 
	size_t dist_tss ; 
	size_t functional ; 	// either INTERGENIC, INTRONIC, or EXONIC

} variant ; 


/* Prints usage */
void usage(const char *pname) ; 

/* Parses the options from command line and returns a struct containing all parameters */
void parse_arguments(parameters* p, char* filename1, char* filename2, int argc, char **argv) ;

/* Prints the parameters to command line */
void print_parameters(parameters* p) ; 

/* Reads in one line of unknown length */
char *inputString(FILE* fp) ; 

/* Returns variant info (when available) */
size_t obtainVariant(FILE *fp, variant *v) ; 

/* Prints the variant to the command line */
void printVariant(variant v) ; 

#endif
