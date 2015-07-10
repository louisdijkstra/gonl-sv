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

__author__ = "Louis Dijkstra"

"""
	The class HumanOmniReader helps with reading the HumanOmni file. 
"""


def strip(string):
	string = string.lstrip()
	return string.rstrip()

class HumanOmniReader:
	
	def __init__(self, location = '/ufs/dijkstra/Projects/SNPs_CNVs/data/HumanOmni25-8v1-1_A-b37.strand'):
		self.file = open(location, 'r') 
		self.positions = [] # list of SNP positions. List is filled when the file is read
	
	def read(self, chromosome):
		"""Reads all the hit SNPs on one chromosome. This must be an integer. 
 		   NOTE: this function must be called first before any other function can be called."""
		self.positions = [] # reset
		chromosome = str(chromosome)
		for line in self.file.readlines():
			raw_data = [value for value in line.split()]
			if strip(raw_data[1]) != chromosome:
				continue
			self.positions.append(int(raw_data[2]))

	def snpPresent(self, position):
		if position in self.positions:
			return True
		return False
			
	def print(self):
		for position in self.positions:
			print('HumanOmni SNP @', position)

	def close(self):
		self.file.close()
