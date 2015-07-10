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



def returnGamma (snp_phase, del_phase, snp_af, del_af):
	g = gamma(snp_phase, del_phase, snp_af, del_af)
	g_n = 0.0 
	if g <= 0:
		return [g, g / snp_af]
	else:
		if del_af <= snp_af:
			return [g, g / (1 - snp_af)]
		else:
			return [g, g / ((snp_af / del_af) - snp_af)] 

def gamma_n (snp_phase, del_phase, snp_af, del_af):
	g = gamma(snp_phase, del_phase, snp_af, del_af)
	if g <= 0:
		return g / snp_af
	else:
		if del_af <= snp_af:
			return g / (1 - snp_af)
		else:
			return g / ((snp_af / del_af) - snp_af) 

def gamma (snp_phase, del_phase, snp_af, del_af):
	p11 = 0.0
	for i in range(len(snp_phase)):
		if snp_phase[i] == 1 and del_phase[i] == 1:
			p11 = p11 + 1.0
	p11 = p11 / len(snp_phase)
	return p11 / del_af - snp_af

