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

def main():
	for chromosome in range(4,23):
		chromosome = str(chromosome)
		print("Starting with chromosome ", chromosome)
		output_file = "results.chr" + chromosome + ".parents.txt"
		command = "python assess_association.py " + chromosome + " > " + output_file
		print('Executing: ', command)
		os.system(command)
		print("Done with chromosome " + chromosome + '. Results are stored in ', output_file)

if __name__ == '__main__':
	sys.exit(main())
