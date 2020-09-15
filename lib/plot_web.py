'''
plot_web : handfull of functions which quickly output the cosmic web features for a given 
snapshot
'''

import pandas as pd 
import upskl_dist as ud 
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os 
import process_segments as ps


def plot_crits(run, runfull, box_side_length, tracer_type, tracer_mass, persistence, snap, ax, subbox = np.array([[0, 75000], [0, 75000], [0, 75000]]), **kwargs):
	'''
	Given a set of parameters defining the chosen disperse run, this plots the set of 
	critical points.

	--- 
	Input : 

	run : str
		Simulation run in format 'TNG100'

	runfull : str
		Simulation run in format 'L75n1820TNG'

	box_side_length : str 
		Simulation box side length in code units. (e.g. TNG100 is 75000)

	tracer_type : str
		Tracer type used to identify cosmic web (e.g. 'stel')

	tracer_mass : str
		Tracer (minimum) mass cutoff for tracers (in format '8-5' i.e. 10**8.5)

	persistence : float
		Persistence ratio used in disperse to consider the robustness of critical 
		point pairs.

	snap : float 
		Snapshot number
	
	ax : plt.axis object
		Axis to plot points onto

	subbox : array 2x3
		Defines limits in x, y and z directions for points to plot.
	
	'''
	 # defining path for critical points file
	crits_filepath = '/Users/chrisduckworth/astronomy/projects/disperse_TNG/upskl_catalogues/' + run + '/' + tracer_type + '_subhalo/crits/'
	crits_filename = run + '_S' + str(snap) + '_M' + tracer_mass + '_' + tracer_type.upper() + '.ascii.NDnet_s' + str(persistence) + '.up.NDskl.BRK.S001.a.crits'

	# loading in ascii file
	crits = ud.read_upskl(crits_filepath + crits_filename)
	# intermediates corresponding to 1-order and 2-order saddle points. 
	minima, saddle_1, saddle_2, nodes = [crits[crits.type == i][['X0','X1','X2']].values for i in np.arange(4)]

	# applying subbox mask to each set of points.
	nmask = nodes[((nodes[:, 0] > subbox[0, 0]) & (nodes[:, 0] <= subbox[0, 1])) & 
				((nodes[:, 1] > subbox[1, 0]) & (nodes[:, 1] <= subbox[1, 1])) & 
				((nodes[:, 2] > subbox[2, 0]) & (nodes[:, 2] <= subbox[2, 1]))]

	mmask = minima[((minima[:, 0] > subbox[0, 0]) & (minima[:, 0] <= subbox[0, 1])) & 
					((minima[:, 1] > subbox[1, 0]) & (minima[:, 1] <= subbox[1, 1])) & 
					((minima[:, 2] > subbox[2, 0]) & (minima[:, 2] <= subbox[2, 1]))]

	s1mask = saddle_1[((saddle_1[:, 0] > subbox[0, 0]) & (saddle_1[:, 0] <= subbox[0, 1])) &
					((saddle_1[:, 1] > subbox[1, 0]) & (saddle_1[:, 1] <= subbox[1, 1])) & 
					((saddle_1[:, 2] > subbox[2, 0]) & (saddle_1[:, 2] <= subbox[2, 1]))]

	s2mask = saddle_2[((saddle_2[:, 0] > subbox[0, 0]) & (saddle_2[:, 0] <= subbox[0, 1])) &
					((saddle_2[:, 1] > subbox[1, 0]) & (saddle_2[:, 1] <= subbox[1, 1])) & 
					((saddle_2[:, 2] > subbox[2, 0]) & (saddle_2[:, 2] <= subbox[2, 1]))]
					  
	# Plotting 
	ax.scatter(nmask[:,0], nmask[:,1], nmask[:,2], c='r', s=10, marker='s', **kwargs)
	ax.scatter(mmask[:,0], mmask[:,1], mmask[:,2], c='dodgerblue', s=10, marker='*', **kwargs)
	ax.scatter(s1mask[:,0], s1mask[:,1], s1mask[:,2], c='g', s=5, marker='^', **kwargs)
	ax.scatter(s2mask[:,0], s2mask[:,1], s2mask[:,2], c='g', s=5, marker='v', **kwargs)
	return 


def plot_segs(run, runfull, box_side_length, tracer_type, tracer_mass, persistence, snap, ax, subbox = np.array([[0, 75000], [0, 75000], [0, 75000]]), **kwargs):
    '''
	Given a set of parameters defining the chosen disperse run, this plots the set of 
	critical points.
    
    --- 
    Input : 
    
    run : str
        Simulation run in format 'TNG100'
    
    runfull : str
        Simulation run in format 'L75n1820TNG'
    
    box_side_length : str 
        Simulation box side length in code units. (e.g. TNG100 is 75000)
    
    tracer_type : str
        Tracer type used to identify cosmic web (e.g. 'stel')
    
    tracer_mass : str
        Tracer (minimum) mass cutoff for tracers (in format '8-5' i.e. 10**8.5)
    
    persistence : float
        Persistence ratio used in disperse to consider the robustness of critical 
        point pairs.
    
    snap : float 
        Snapshot number
        
    ax : plt.axis object
    	Axis to plot points onto
    
    subbox : array 2x3
    	Defines limits in x, y and z directions for points to plot.
    	
    '''
    # defining path for filament segment file
    segs_filepath = '/Users/chrisduckworth/astronomy/projects/disperse_TNG/upskl_catalogues/' + run + '/' + tracer_type + '_subhalo/segs/'
    segs_filename = run + '_S' + str(snap) + '_M' + tracer_mass + '_' + tracer_type.upper() + '.ascii.NDnet_s' + str(persistence) + '.up.NDskl.BRK.S001.a.segs'

    # extracting segments from segs tabledata.
    segs = ud.read_upskl(segs_filepath + segs_filename)
    # Extracting start and end of segments. Calculating middpoints within periodic box.
    U = segs[['U0','U1','U2']].values
    V = segs[['V0','V1','V2']].values
    # finding midpoints since this is quicker to plot than individual lines.
    segs_mid = ud.seg_midpoint(U, V, box_side_length, periodic=False)
    
    segmask = segs_mid[((segs_mid[:, 0] > subbox[0, 0]) & (segs_mid[:, 0] <= subbox[0, 1])) & 
    					((segs_mid[:, 1] > subbox[1, 0]) & (segs_mid[:, 1] <= subbox[1, 1])) & 
    					((segs_mid[:, 2] > subbox[2, 0]) & (segs_mid[:, 2] <= subbox[2, 1]))]
                   
    ax.scatter(segmask[:,0], segmask[:,1], segmask[:,2], c='k', s=2, marker='o', **kwargs) 
    return 