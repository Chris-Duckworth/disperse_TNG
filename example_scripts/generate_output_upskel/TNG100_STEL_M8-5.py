'''
Generates FOF and subhalo distances for TNG100_STEL_M8-5 

run = 'TNG100'
runfull = 'L75n1820TNG'
box_side_length = 75000
tracer_type = 'stel'
tracer_mass = '8-5'
persistence = 4
'''

# -----------------------------------------------------------------------------
import numpy as np
import pandas as pd
import upskl_dist as ud
import os
from scipy import spatial
import h5py

# -----------------------------------------------------------------------------
# Box parameters

run = 'TNG100'
runfull = 'L75n1820TNG'
box_side_length = 75000
snapshots = np.array([33, 40, 50, 67, 78, 84, 91, 99]) #np.arange(30, 100)
tracer_type = 'stel'
tracer_mass = '8-5'
persistence = 4

# -----------------------------------------------------------------------------
# For all snapshots.

for snap in snapshots:
	ud.create_cw_file(run, runfull, box_side_length, tracer_type, tracer_mass, persistence, snap)