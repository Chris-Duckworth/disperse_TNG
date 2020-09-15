'''
upskl_dist.py

Given a defined skeleton (critical points [.crits] and segments [.segs]) and a set of object positions (e.g. galaxies), this script returns the distance to the nearest cosmic web feature of each type (minima, saddle, node, filament) in the up-skeleton.

Chris J Duckworth cduckastro@gmail.com
'''

import numpy as np
import itertools
import pandas as pd
from scipy import spatial
import h5py
import os 


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


def create_cw_file(run, runfull, box_side_length, tracer_type, tracer_mass, persistence, snap):
    '''
    Given a set of params defining where to look-up the raw disperse files, this creates and saves a hdf5 file 
    containing the distances to cosmic web features (as defined by disperse) for all fofs and subhaloes in the 
    snapshot.
    
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
    '''
    
    # defining path for critical points file
    crits_filepath = '/Users/chrisduckworth/astronomy/projects/disperse_TNG/upskl_catalogues/' + run + '/' + tracer_type + '_subhalo/crits/' + 'M' + tracer_mass + '/'
    crits_filename = run + '_S' + str(snap) + '_M' + tracer_mass + '_' + tracer_type.upper() + '.ascii.NDnet_s' + str(persistence) + '.up.NDskl.BRK.S001.a.crits'

    # loading in ascii file
    crits = read_upskl(crits_filepath + crits_filename)
    # intermediates corresponding to 1-order and 2-order saddle points. 
    minima, saddle_1, saddle_2, nodes = [crits[crits.type == i][['X0','X1','X2']].values for i in np.arange(4)]
    
    # defining path for filament segment file
    segs_filepath = '/Users/chrisduckworth/astronomy/projects/disperse_TNG/upskl_catalogues/' + run + '/' + tracer_type + '_subhalo/segs/' + 'M' + tracer_mass + '/'
    segs_filename = run + '_S' + str(snap) + '_M' + tracer_mass + '_' + tracer_type.upper() + '.ascii.NDnet_s' + str(persistence) + '.up.NDskl.BRK.S001.a.segs'

    # extracting segments from segs tabledata.
    segs = read_upskl(segs_filepath + segs_filename)
    # Extracting start and end of segments. Calculating middpoints within periodic box.
    U = segs[['U0','U1','U2']].values
    V = segs[['V0','V1','V2']].values
    segs_mid = seg_midpoint(U, V, box_side_length, periodic=False)
    
    # Creating KDtrees for each set of critical points and segment midpoints.
    print('Generating KDtree for CW positions.')
    d_node_tree = spatial.cKDTree(nodes)
    d_saddle_1_tree = spatial.cKDTree(saddle_1)
    d_saddle_2_tree = spatial.cKDTree(saddle_2)
    d_minima_tree = spatial.cKDTree(minima)
    d_skel_tree = spatial.cKDTree(segs_mid)
    
    # loading in all FoFs and subhaloes
    lookup_path = '/Users/chrisduckworth/astronomy/projects/disperse_TNG/tracer_catalogues/' + run + '/general_halo_subhalo/'
    fof_name = runfull + '-SNAP' + str(snap) + '-HALOS.hdf5'
    subhalo_name = runfull + '-SNAP' + str(snap) + '-SUBHALOS.hdf5'

    fof = h5py.File(lookup_path + fof_name, 'r')
    subhalo = h5py.File(lookup_path + subhalo_name, 'r')
    print(fof_name+' loaded.')
    print(subhalo_name+' loaded.')
    
    fof_pos = fof['pos'][()]
    subhalo_pos = subhalo['pos'][()]
    
    # finding cw distances for fofs.
    d_node_fof = d_node_tree.query(fof_pos)[0]
    d_saddle_1_fof = d_saddle_1_tree.query(fof_pos)[0]
    d_saddle_2_fof = d_saddle_2_tree.query(fof_pos)[0]
    d_minima_fof = d_minima_tree.query(fof_pos)[0]
    d_skel_fof = d_skel_tree.query(fof_pos)[0]
    
    # finding cw distances for subhaloes.
    d_node_subhalo = d_node_tree.query(subhalo_pos)[0]
    d_saddle_1_subhalo = d_saddle_1_tree.query(subhalo_pos)[0]
    d_saddle_2_subhalo = d_saddle_2_tree.query(subhalo_pos)[0]
    d_minima_subhalo = d_minima_tree.query(subhalo_pos)[0]
    d_skel_subhalo = d_skel_tree.query(subhalo_pos)[0]
    
    # defining outpath and file name structure
    outpath = '/Users/chrisduckworth/astronomy/projects/disperse_TNG/output_upskl/' + run + '/' + tracer_type + '_subhalo/'
    outfile = run + '_S' + str(snap) + '_M' + tracer_mass + '_' + tracer_type.upper()

    # creating output fof file.
    print('Creating output : ' + outpath + 'fof_' + outfile + '.hdf5')
    fof_output = h5py.File(outpath + 'fof_' + outfile + '.hdf5', 'w')

    # adding columns to fof file.
    fof_output.create_dataset('group_ID', data=np.arange(fof_pos.shape[0]))
    #fof_output.create_dataset('pos', data=fof_pos)
    fof_output.create_dataset('d_node', data=d_node_fof)
    fof_output.create_dataset('d_saddle_1', data=d_saddle_1_fof)
    fof_output.create_dataset('d_saddle_2', data=d_saddle_2_fof)
    fof_output.create_dataset('d_minima', data=d_minima_fof)
    fof_output.create_dataset('d_skel', data=d_skel_fof)
    fof_output.close()

    # creating output subhalo file
    print('Creating output : ' + outpath + 'subhalo_' + outfile + '.hdf5')
    subhalo_output = h5py.File(outpath + 'subhalo_' + outfile + '.hdf5', 'w')

    # adding columns to subhalo file
    subhalo_output.create_dataset('subhalo_ID', data=np.arange(subhalo_pos.shape[0]))
    #subhalo_output.create_dataset('pos', data=subhalo_pos)
    subhalo_output.create_dataset('d_node', data=d_node_subhalo)
    subhalo_output.create_dataset('d_saddle_1', data=d_saddle_1_subhalo)
    subhalo_output.create_dataset('d_saddle_2', data=d_saddle_2_subhalo)
    subhalo_output.create_dataset('d_minima', data=d_minima_subhalo)
    subhalo_output.create_dataset('d_skel', data=d_skel_subhalo)
    subhalo_output.close()
    return 







