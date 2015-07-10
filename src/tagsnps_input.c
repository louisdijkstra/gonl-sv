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

#include "tagsnps_input.h"

void usage(const char *pname) {
	printf(	"Usage: %s [OPTION] <deletions.haplotypes> <snps.haplotypes>\n"
		"\n"
		" -p\tNUM\tP-value threshold.\n"
		" -r\tNUM\tR^2 threshold.\n"
		" -s\tNUM\tRange in which to search for appropriate tag SNPs. \n"
		" -v\t\tVerbose.\n"
		" -h\t\tPrint this help.\n\n", 
		pname) ; 
	exit(0) ; 
}

void print_parameters(parameters* p) {
	printf(	"===PARAMETER SETTINGS===\n\n"
		"p-value thres\t%lf\n"
		"R^2 thres\t\t%lf\n"
		"search range\t%zd\n"
		"verbose\t\t%zd\n\n",
		p->p_threshold,  
		p->r2_threshold, 
		p->search_range,
		p->verbose
		) ; 
}

void parse_arguments(parameters* p, char* input_filename1, char* input_filename2, int argc, char **argv) {
	
	/* set defaults */
	p->p_threshold = 0.05 ;   
	p->r2_threshold = 0.8 ;  
	p->search_range = 1000000 ; 
	p->verbose	= 0 ; 
	
	size_t ch ; 
	while ((ch = getopt(argc, argv, "p:r:s:vh")) != -1)
    	{
       		switch(ch) {
			case 'p': p->p_threshold = strtod(optarg, 0); break ; 
			case 'r': p->r2_threshold = strtod(optarg, 0); break ; 
			case 's': p->search_range = strtol(optarg, 0, 10); break ; 
			case 'v': p->verbose = 1; break ; 
        		case 'h': default: usage(argv[0]);
        	}
    	}

	if (argv[optind] == NULL || argv[optind+1] == NULL) {
		usage(argv[0]) ; 
		exit(EXIT_SUCCESS) ; 
	}
	
	// check whether parameter settings are valid 
	size_t error_occurred = 0 ; 

	if (p->p_threshold < 0 || p->p_threshold > 1.0) {
		printf("ERROR: Invalid argument. p-value threshold (-p) must lie in the unit interval.\n") ;
		error_occurred = 1 ;  
	}

	if (p->r2_threshold < 0 || p->r2_threshold > 1.0) {
		printf("ERROR: Invalid argument. r2 threshold (-r) must lie in the unit interval.\n") ;
		error_occurred = 1 ;  
	}

	if (p->search_range < 0) {
		printf("ERROR: Invalid argument. search_range (-s) must be nonnegative.\n") ;
		error_occurred = 1 ;  
	}

	if (error_occurred == 1) {
		exit(EXIT_FAILURE) ; 
	}

	sprintf(input_filename1, argv[optind]) ;
	sprintf(input_filename2, argv[optind+1]) ;
}

char *inputString(FILE* fp) {
	/* Reads in one line of unknown length */
	size_t buffer_length = BUFSIZ ;
	char *buffer = malloc(buffer_length * sizeof(char)) ; 
	if (buffer == NULL) {
		printf("Insufficient memory to allocate the buffer.\n") ;
		exit(EXIT_FAILURE) ; 
	}

	int ch ; 
	size_t i = 0; 
	while (EOF != (ch=fgetc(fp)) && ch != '\n') {
		buffer[i++] = ch ; 
		if (i == buffer_length) {
			buffer = realloc(buffer, sizeof(char) * (buffer_length += BUFSIZ)) ; 
			if (buffer == NULL) {
				printf("Insufficient memory to reallocate the buffer.\n") ;
				exit(EXIT_FAILURE) ; 
			}
		}
	}
	buffer[i++] = '\0' ; 
	return realloc(buffer, sizeof(char) * i); 
}

size_t obtainVariant(FILE *fp, variant *v) {
	/* Returns a variant and stores it in v. In case there is no variant no more, 
	 * the function returns END_OF_FILE_REACHED. Otherwise, it outputs 
	 * VALID_VARIANT. 
	 */
	char ch ;
	size_t i ; 
	char length ; 
	if(fscanf(fp, "%zd %c %zd%c", &(v->autosome), &(v->type), &(v->position), &ch) != 4) {
		return END_OF_FILE_REACHED ; 
	}
	if (v->type == '+' || v->type == '-') {
		fscanf(fp, "%zd%c", &(v->length), &ch) ; 
	} else {
		fscanf(fp, "%c%c", &length, &ch) ; 
	}

	fscanf(fp, "%c %c %zd %lf%c", &(v->ref), &(v->alt), &(v->n), &(v->af), &ch) ; 
	
	char * buffer = malloc(v->n * sizeof(char)) ; 
	v->haplotypes = malloc(v->n * sizeof(size_t)) ; 
	if (v->haplotypes == NULL || buffer == NULL) {
		printf("ERROR: insufficient memory to allocate haplotype data.\n") ; 
		exit(EXIT_FAILURE) ; 	
	}
	fscanf(fp, "%s%c", buffer, &ch) ; 
	for (i = 0; i < v->n; i ++) {
		v->haplotypes[i] = buffer[i] - '0' ;
	}

	if (ch == '\n') {
		return VALID_VARIANT ; 
	}

	// read till the end of the line
	do {
		ch = fgetc(fp);
	} while (ch != '\n') ; 

	return VALID_VARIANT ; 
}

void printVariant (variant v) {
	/* Prints the data to the command line */
	size_t i ; 
	printf("autosome:\t%zd\n", v.autosome) ; 
	printf("position:\t%zd\n", v.position) ; 
	printf("type:\t\t%c\n", v.type) ; 
	if (v.type == '-' || v.type == '+') {
		printf("length:\t%zd\n", v.length) ; 
	}
	printf("# observations: %zd\tfrequency: %f\n", v.n, v.af) ; 
	for (i = 0; i < v.n; i ++) {
		printf("%zd", v.haplotypes[i]) ; 
	}
	printf("\n\n") ; 
}

