'''
select_subhalos

library of functions for selecting subhalos (for purpose of cw reconstruction) and padding 
set of positions in periodic cube (to counteract disperse's problems with boundary conds.)

Chris J Duckworth cduckastro@gmail.com
'''

import numpy as np 
import groupcat as gc


def return_stel_tracers(basepath, snapnum, blen, min_mass=10**8.5):
	'''
	returns stellar tracer positions for all objects with a defined minimum mass at a 
	defined snapshot. basepath should direct to /output directory of chosen simulation 
	and blen should be the corresponding box side length.
	'''
	# Subhalo fields
	sfields = ['SubhaloMassType','SubhaloPos']
	subhalos = gc.loadSubhalos(basepath, snapnum, fields=sfields)
	
	# selecting based on total mass within the subhalo.
	return subhalos['SubhaloPos'][()][subhalos['SubhaloMassType'][()][:,4]*10**10 > min_mass]


def box_extend(blen, pos, frac=0.1):
	'''
	Given a periodic cube of equal side length blen, this script extends (or pads) each
	side by an additional fraction set by frac. This is simply repeating all 
	positions (defined by pos) on each side.
	
	This is required to ensure disperse is doing its job at the edges of sims.
	'''
	
	# wrapping x, y and then z directions in order. wrapping is performed on updated set 
	# of positions at each step to ensure there is enough padding at the corners.
	
	# creating copy of positions to update.
	pos_extend = pos

	# looping over dimensions in order.
	for ind in np.arange(3):
		# slicing each side.
		low_slice = pos_extend[pos_extend[:, ind] < blen * frac] 
		low_slice[:, ind] += blen 
	
		high_slice = pos_extend[pos_extend[:, ind] > blen * (1 - frac)]
		high_slice[:, ind] -= blen
		
		# adding new points.
		pos_extend = np.append(pos_extend, low_slice, axis=0)
		pos_extend = np.append(pos_extend, high_slice, axis=0)
	
	return pos_extend

