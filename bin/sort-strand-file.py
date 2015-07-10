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
import math
import vcf

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))[:-3] + 'python')
from Variant import *

__author__ = "Louis Dijkstra"

usage = """%prog <strand-file>

Reads in a sorts a strand file (e.g., human omni data)
"""

def main():
	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()

	if (len(args)!=1):
		parser.print_help()
		return 1

	strand_file = open(args[0], 'r')

	snp_positions = [[] for i in range(1,23)] # one for every autosome
	
	for line in strand_file:
		values = line.split('\t')
		autosome = values[1].strip()
		if autosome == 'MT' or autosome == 'Y' or autosome == 'X':
			continue	
		snp_positions[int(autosome) - 1].append(int(values[2])) 

	# sort array of positions per autosome
	for autosome in range(22):
		snp_positions[autosome].sort()
	
	# print results to screen
	for autosome in range(22):
		for i in range(len(snp_positions[autosome])):
			print("%d %d"%(autosome+1, snp_positions[autosome][i]))

if __name__ == '__main__':
	sys.exit(main())

