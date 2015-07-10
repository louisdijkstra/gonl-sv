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
import sys
import os 

__author__ = "Louis Dijkstra"

"""
	Tags all deletions per chromosome
"""

def main():
	# list of autosomes to be processed 
	list_of_autosomes = [11,10,9,8,7,6,5,4,3,1] # TODO change back if necessary
	for autosome in list_of_autosomes:
		autosome = str(autosome)
		file_location = 'Results/tags_chr' + autosome + '.txt'
		command = 'python tag_deletions.py -c ' + autosome + ' > ' + file_location
		print('Executing: ', command)
		os.system(command)

if __name__ == '__main__':
	sys.exit(main())
