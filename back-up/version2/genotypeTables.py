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

import math
from scipy.stats.stats import pearsonr

def returnN (g):
	sum(g)

def rowTotals (g):
	rt = [0,0,0]
	rt[0] = g[0] + g[1] + g[2]
	rt[1] = g[3] + g[4] + g[5]
	rt[2] = g[6] + g[7] + g[8]
	return rt 

def columnTotals (g):
	ct = [0,0,0]
	ct[0] = g[0] + g[3] + g[6] 
	ct[1] = g[1] + g[4] + g[7]
	ct[2] = g[2] + g[5] + g[8]  
	return ct 


def boundR (g):
	H = float(2 * sum(g)) # total number of observed haplotypes
	rt = rowTotals(g)
	ct = columnTotals(g)
	h_A = 2 * ct[0] + ct[1]
	h_B = 2 * rt[0] + rt[1]
	p_A = h_A / H 
	p_B = h_B / H
	norm = math.sqrt(p_A * (1 - p_A) * p_B * (1 - p_B)) 
	print norm 
	return [D(g,0.0) / norm, D(g,1.0) / norm]	

def boundD (g):
	return [D(g, 0.0), D(g,1.0)]


def D(g, pi):
	H = float(2 * sum(g)) # total number of observed haplotypes
	rt = rowTotals(g)
	ct = columnTotals(g)
	h_A = 2 * ct[0] + ct[1]
	h_B = 2 * rt[0] + rt[1]
	D = (2 * g[0] + g[1] + g[3] + g[4] * float(pi)) / H - (h_A * h_B) / H**2
	return D

def print3x3Table (g):

	rt = rowTotals(g)
	ct = columnTotals(g)	

	print '\t11\t10/01\t00\t|\ttotal' 
	print '-----------------------------------------------'
	print '11\t', g[0], '\t', g[1], '\t', g[2], '\t|\t', rt[0]
	print '10/01\t', g[3], '\t', g[4], '\t', g[5], '\t|\t', rt[1]
	print '00\t', g[6], '\t', g[7], '\t', g[8], '\t|\t', rt[2] 
	print '-----------------------------------------------'	
	print '\t', ct[0], '\t', ct[1], '\t', ct[2], '\t|\t', sum(g) 


def r3x3Table (x,y):
	s = []
	t = []
	for i in range(0, len(x), 2):
		s.append(abs(x[i] + x[i+1] - 2))
		t.append(abs(y[i] + y[i+1] - 2))
	return pearsonr(s,t)
	

def reconstruct3x3Table (x,y):
	g = [0] * 9 
	for i in range(0, len(x), 2):
		s = abs(x[i] + x[i+1] - 2)
		t = abs(y[i] + y[i+1] - 2)
		g[t * 3 + s] = g[t * 3 + s] + 1
	return g 
 



