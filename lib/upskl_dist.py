'''
upskl_dist.py

Given a defined skeleton (critical points [.crits] and segments [.segs]) and a set of object positions (e.g. galaxies), this script returns the distance to the nearest cosmic web feature of each type (minima, saddle, node, filament) in the up-skeleton.

Chris J Duckworth cduckastro@gmail.com
'''

import numpy as np
import itertools
import pandas as pd
from scipy import spatial


def read_upskl(path_to_file):
	'''
	reads in .crits or .seg files directly from disperse output and returns a pandas 
	dataframe object
	'''
	with open(path_to_file) as handle:
        # Reading in all lines starting with # and selecting the second row which has column names.
		_, names, *_ = itertools.takewhile(lambda line: line.startswith('#'), handle)
		names = names[1:].split()
	return pd.read_csv(path_to_file, header=0, names=names, sep=' ', comment='#')


def nearest_neighbour(lookup_positions, upskl_positions):
	'''
	constructs KDtree for set of cosmic web positions (i.e. specific set of critical points
	or segment midpoints and then returns the distance to the nearest for each of the 
	supplied look-up positions.
	
	DISCLAIMER: this is not necessarily faster than brute force for small numbers of 
	query positions.
	'''
	tree = spatial.cKDTree(upskl_positions)
	return tree.query(lookup_positions)


def seg_midpoint(U, V, periodic=False, blen=75000):
	'''
	provides midpoints for all segments defined by U and V.
	
	to get around disperse failing to reconnect periodic box sides when sampling is low, 
	the box has been extended in each dimension. therefore segments will be defined outside
	of the periodic cube. this is fine if you are only interested in distances to various 
	cosmic web features, however, continue with trepidation if you want to investigate 
	properties of the cosmic web itself since some will be more or less repeated.
	'''
	
	if periodic == False:
		print('assuming box has been extended and hence no-longer has periodic boundaries!')
		return (U + V)/2
		
	elif periodic == True:
		print('assuming periodic cube - midpoints will be returned in range defined by box length')
		# Using the periodic distance calculator to identify segments going across the boundary.
		del1 = np.abs(U - V)
		del2 = blen - np.abs(U - V)
		segs_mid = np.zeros(del1.shape)
		segs_mid[del1 <= del2] = (U[del1 <= del2] + V[del1 <= del2]) / 2
		# Selecting all segments which cross the boundary. The midpoint is estimated by adding the box length in the average.
		segs_mid[del1 > del2] = (U[del1 > del2] + V[del1 > del2] + blen) / 2
		# Sometimes the midpoint will be above the boundary edge. Taking off box length to fix if this happens.
		segs_mid[segs_mid > blen] = segs_mid[segs_mid > blen] - blen
		return segs_mid
		
	else :
		assert False, "Periodic must be true or false"

