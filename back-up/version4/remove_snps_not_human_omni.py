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

from __future__ import print_function, division
from optparse import OptionParser
import os
import sys

from SNPDeletionPair import * 
from HumanOmniReader import * 

__author__ = "Louis Dijkstra"

usage = """%prog <results_file> <chromosome>

	<results_file> 	- file with results as produced by assess_association.py. 
	<chromosome> 	- the chromosome in question (integer)

Removes every SNP-deletion pair from the file for which the SNP is not present on the HumanOmni array. 
"""

def main():
	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()

	if (len(args)!=2):
		parser.print_help()
		return 1

	results_file	 = open(args[0], 'r')
	chromosome 	 = int(args[1])

	human_omni_reader = HumanOmniReader()
	human_omni_reader.read(chromosome)

	for line in results_file:
		pair = SNPDeletionPair(line)
		if not human_omni_reader.snpPresent(pair.snp.position):
			continue
		print(line, end = '')

	results_file.close()

if __name__ == '__main__':
	sys.exit(main())
