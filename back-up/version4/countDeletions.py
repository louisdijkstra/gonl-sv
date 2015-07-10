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
import scipy as sp
import scipy.stats
from scipy.stats.stats import pearsonr

from DeletionReader import *

__author__ = "Louis Dijkstra"

usage = """%prog

Counts the number of GoNL deletions with a minor allele frequency above the 4%
"""

def main():

	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()

	vcf_reader = vcf.Reader(open("/data1/gonl/new-deletion-phasing/deletions-20-10000-phased-20140925.vcf"))

	n_del = 0 
	for vcf_record in vcf_reader:
		deletion = Deletion(vcf_record)
		if deletion.major_af > .96:
			continue
		n_del += 1
		if n_del % 10 == 0:
			print(n_del, ' deletions processed so far.')
	print("Total number of deletions: ", n_del)
	

if __name__ == '__main__':
	sys.exit(main())
