'''
general_halo_subhalo_files

Simply grabbing halo and subhalo files with positions.
'''

import numpy as np 
import groupcat as gc
import h5py 

# ---------------------------------------------------------------------------------------
# Defining for fiducial TNG100.

TNGrun, blen = 'L75n1820TNG', 75000  # comoving kpc/h
basepath = '/virgo/simulations/IllustrisTNG/'+TNGrun+'/output' 
outpath = '/u/cduckwor/python/disperse_release/tracer_catalogues/TNG100/general_halo_subhalo/'

# ---------------------------------------------------------------------------------------
# Grabbing all haloes and subhaloes.

snaps = np.arange(100) #np.array([33, 50, 67, 72, 78, 84, 91, 99]) 

for snapnum in snaps:
	# loading group positions
	gfields = ['GroupPos']
	halos = gc.loadHalos(basepath, snapnum, fields=gfields)
	halo_name = TNGrun + '-SNAP' + str(snapnum) + '-HALOS.hdf5'

	# Saving halo file 
	halo_hf = h5py.File(outpath + halo_name, 'w')
	halo_hf.create_dataset('pos', data = halos )
	halo_hf.close()


	# loading halo positions
	sfields = ['SubhaloPos']
	subhalos = gc.loadSubhalos(basepath, snapnum, fields=sfields)
	subhalo_name = TNGrun + '-SNAP' + str(snapnum) + '-SUBHALOS.hdf5'

	# Saving subhalo file
	subhalo_hf = h5py.File(outpath + subhalo_name, 'w')
	subhalo_hf.create_dataset('pos', data = subhalos )
	subhalo_hf.close()

# ---------------------------------------------------------------------------------------
# Defining for fiducial TNG300.

TNGrun, blen = 'L205n2500TNG', 205000  # comoving kpc/h
basepath = '/virgo/simulations/IllustrisTNG/'+TNGrun+'/output' 
outpath = '/u/cduckwor/python/disperse_release/tracer_catalogues/TNG300/general_halo_subhalo/'

# ---------------------------------------------------------------------------------------
# Grabbing all haloes and subhaloes.

snaps = np.arange(100) #np.array([33, 50, 67, 72, 78, 84, 91, 99]) 

for snapnum in snaps:
	# loading group positions
	gfields = ['GroupPos']
	halos = gc.loadHalos(basepath, snapnum, fields=gfields)
	halo_name = TNGrun + '-SNAP' + str(snapnum) + '-HALOS.hdf5'

	# Saving halo file 
	halo_hf = h5py.File(outpath + halo_name, 'w')
	halo_hf.create_dataset('pos', data = halos )
	halo_hf.close()


	# loading halo positions
	sfields = ['SubhaloPos']
	subhalos = gc.loadSubhalos(basepath, snapnum, fields=sfields)
	subhalo_name = TNGrun + '-SNAP' + str(snapnum) + '-SUBHALOS.hdf5'

	# Saving subhalo file
	subhalo_hf = h5py.File(outpath + subhalo_name, 'w')
	subhalo_hf.create_dataset('pos', data = subhalos )
	subhalo_hf.close()

# ---------------------------------------------------------------------------------------
