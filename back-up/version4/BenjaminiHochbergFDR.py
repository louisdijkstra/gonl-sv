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
	Performs Benjamini and Hochbergs false discovery rate control procedure on 
	a list of p-values and a list of associated objects. It returns the list
	of objects that are considered significant. These objects can be realizations 
	of a certain class, for example.
"""

def benjamini_hochberg (q, p_values, objects):
	"""Returns the list of objects that are deemed significant."""
	significant_objects = [] 	

	# sort the p_values and the objects list accordingly. 
	p_values, objects = (list(x) for x in zip(*sorted(zip(p_values, objects))))
	
	number_of_tests = len(p_values)
	for k in range(1, number_of_tests + 1):
		if p_values[k-1] <= float(k) / float(number_of_tests) * q:
			significant_objects.append(objects[k-1])
		else:
			break 
	return significant_objects 
	
		


