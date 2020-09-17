'''
DM_subhalo_sel.py 

This script selects a DM set of tracers (i.e. DM subhaloes with a minimum 
mass of 10**8.5 Msol) for comparison to observations.
'''

import numpy as np 
import h5py 
import groupcat as gc
#import coordinate_transforms as ct
import matplotlib.pyplot as plt
import select_subhalos as sel
import pandas as pd

# ---------------------------------------------------------------------------------------
# Defining for fiducial TNG100.

TNGrun, blen = 'L75n1820TNG', 75000  # comoving kpc/h
basepath = '/virgo/simulations/IllustrisTNG/'+TNGrun+'/output' 
outpath = '/u/cduckwor/python/disperse_release/tracer_catalogues/TNG100/DM_subhalo/'

# ---------------------------------------------------------------------------------------
# running for a defined set of snapshots.

snaps = np.array([40, 59]) # np.array([33, 40, 50, 59, 67, 72, 78, 84, 91, 99]) # np.arange(100)

for ind in snaps:
	# wrapping 10% of box edges in each direction.
	pos = sel.return_DM_tracers(basepath, ind, blen, min_mass=10**8)
	pos_wrap = sel.box_extend(blen, pos, frac=0.1)
	
	print('Computing for snap-'+str(ind)+'. There are '+str(pos.shape[0])+
		  ' original galaxies and '+str(pos_wrap.shape[0])+' galaxies after box extension.')
	
	# Exporting in disperse format.
	tab = pd.DataFrame(pos_wrap, columns=['#px', 'py', 'pz'])
	outfile = 'TNG100_S'+str(ind)+'_M8_DM.ascii'
	tab.to_csv(outpath+outfile, index=False)

# ---------------------------------------------------------------------------------------
# Also running for TNG300.

TNGrun, blen = 'L205n2500TNG', 205000  # comoving kpc/h
basepath = '/virgo/simulations/IllustrisTNG/'+TNGrun+'/output' 
outpath = '/u/cduckwor/python/disperse_release/tracer_catalogues/TNG300/DM_subhalo/'

# ---------------------------------------------------------------------------------------
# running for a defined set of snapshots.

snaps = np.array([40, 59]) # np.array([33, 40, 50, 59, 67, 72, 78, 84, 91, 99]) # np.arange(100)

for ind in snaps:
	# wrapping 10% of box edges in each direction.
	pos = sel.return_DM_tracers(basepath, ind, blen, min_mass=10**8)
	pos_wrap = sel.box_extend(blen, pos, frac=0.1)
	
	print('Computing for snap-'+str(ind)+'. There are '+str(pos.shape[0])+
		  ' original galaxies and '+str(pos_wrap.shape[0])+' galaxies after box extension.')
	
	# Exporting in disperse format.
	tab = pd.DataFrame(pos_wrap, columns=['#px', 'py', 'pz'])
	outfile = 'TNG300_S'+str(ind)+'_M8_DM.ascii'
	tab.to_csv(outpath+outfile, index=False)

# ---------------------------------------------------------------------------------------
