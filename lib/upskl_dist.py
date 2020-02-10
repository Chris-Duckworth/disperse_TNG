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


def seg_midpoint(U, V, box_side_length, periodic=False):
	'''
	Returns midpoints for all segments defined by U and V (in DisPerSE). 
	
	To get around disperse failing to reconnect periodic box sides when sampling is low, 
	the box has been extended in each dimension. Therefore segments will be defined outside
	of the periodic cube. 
	
	process_segments reconnects segments and re-enforces periodic boundaries, in that case 
	use periodic=True which calculates midpoints taking into account periodic boundaries.
	
	'''
	
	if periodic == False:
		print('assuming box has been extended and hence no-longer has periodic boundaries!')
		return (U + V)/2
		
	elif periodic == True:
		print('assuming periodic cube - midpoints will be returned in range defined by box length')
		# Using the periodic distance calculator to identify segments going across the boundary.
		del1 = np.abs(U - V)
		del2 = box_side_length - np.abs(U - V)
		segs_mid = np.zeros(del1.shape)
		segs_mid[del1 <= del2] = (U[del1 <= del2] + V[del1 <= del2]) / 2
		# Selecting all segments which cross the boundary. The midpoint is estimated by adding the box length in the average.
		segs_mid[del1 > del2] = (U[del1 > del2] + V[del1 > del2] + box_side_length) / 2
		# Sometimes the midpoint will be above the boundary edge. Taking off box length to fix if this happens.
		segs_mid[segs_mid > box_side_length] = segs_mid[segs_mid > box_side_length] - box_side_length
		return segs_mid
		
	else :
		assert False, "Periodic must be true or false"

