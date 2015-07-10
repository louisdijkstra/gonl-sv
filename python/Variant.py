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
import vcf

__author__ = "Louis Dijkstra"

def returnAutosome (autosome):
	"""Returns an integer denoting the chromosome."""
	if len(autosome) > 3 and autosome[:3] == 'chr':
		autosome = autosome[3:]
	if autosome == 'X' or autosome == 'x':
		return 23
	if autosome == 'Y' or autosome == 'y':
		return 24	
	return int(autosome)

def returnIndelLength (vcf_record):
	"""Returns variation length of an indel given its vcf record"""
	if 'SVLEN' in vcf_record.INFO:
		if isinstance(vcf_record.INFO['SVLEN'], list):
			return abs(vcf_record.INFO['SVLEN'][0]) 
		else:
			return abs(vcf_record.INFO['SVLEN']) 
	else:
		return abs(len(vcf_record.REF) + 1 - len(vcf_record.ALT[0]))

def isIndel (vcf_record):
	"""Determines whether a vcf record represents an indel or not. WARNING: only considers the first alternative (ALT[0]); others are neglected."""
	return (isDeletion(vcf_record) or isInsertion(vcf_record))

def isDeletion (vcf_record):
	"""Determines whether a vcf record represents a deletion or not. NOTE: only considers first alternative (ALT[0]); others are neglected"""
	if 'SVTYPE' in vcf_record.INFO: # checks whether the key SVTYPE is present in INFO
		if isinstance(vcf_record.INFO['SVTYPE'], list):
			if vcf_record.INFO['SVTYPE'][0] == 'DEL':	
				return True
		else:
			if vcf_record.INFO['SVTYPE'] == 'DEL':	
				return True
	else: # makes statement on the basis of length of ALT and REF 
		if len(vcf_record.REF) > 1 and len(vcf_record.ALT[0]) == 1:
			return True
	return False

def isInsertion (vcf_record):
	"""Determines whether a vcf record represents an insertion or not. NOTE: only considers first alternative (ALT[0]); others are neglected"""
	if 'SVTYPE' in vcf_record.INFO: # checks whether the key SVTYPE is present in INFO
		if isinstance(vcf_record.INFO['SVTYPE'], list):
			if vcf_record.INFO['SVTYPE'][0] == 'INS':	
				return True
		else:
			if vcf_record.INFO['SVTYPE'] == 'INS':	
				return True
	else: # makes statement on the basis of length of ALT and REF 
		if len(vcf_record.REF) == 1 and len(vcf_record.ALT[0]) > 1:
			return True
	return False

def isSNP (vcf_record):
	"""Determines whether a vcf record represents a SNP or not. NOTE: only considers first alternative (ALT[0]); others are neglected"""
	if len(vcf_record.REF) == 1 and len(vcf_record.ALT[0]) == 1:
		return True
	return False
