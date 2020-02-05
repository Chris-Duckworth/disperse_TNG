'''
find_upskl_distances

Example script which looks up and finds nearest critical points (nodes, saddles and 
minima) and nearest segment.

assumes you have disperse_TNG/lib/ in your pythonpath.
'''

import numpy as np
import pandas as pd
import upskl_dist as ud
import os 

# -----------------------------------------------------------------------------
# defining overall filepath for upskl files (crits and segs).
filepath = os.environ['DISPERSE']+'/disperse_TNG/upskl_catalogues/TNG100/stel_subhalo/'
cname = 'TNG100_S99_M8-5_STEL.ascii.NDnet_s4.up.NDskl.BRK.S001.a.crits' # snapshot of choice.
sname = 'TNG100_S99_M8-5_STEL.ascii.NDnet_s4.up.NDskl.BRK.S001.a.segs' # snapshot of choice.

# defining random 3d positions.
pos = np.random.rand(100,3) * 75000

# -----------------------------------------------------------------------------
# extracting nodes, saddles & minima from crits tabledata. 
crits = ud.read_upskl(filepath+'crits/'+cname)
# intermediates corresponding to 1-order and 2-order saddle points. 
minima, saddle_1, saddle_2, nodes = [crits[crits.type == i][['X0','X1','X2']].values for i in np.arange(4)]

# Creating KDtrees for each set of critical points & returning distance to nearest neighbour.
d_node = ud.nearest_neighbour(pos, nodes)[0]
d_saddle_1 = ud.nearest_neighbour(pos, saddle_1)[0]
d_saddle_2 = ud.nearest_neighbour(pos, saddle_2)[0]
d_minima = ud.nearest_neighbour(pos, minima)[0]

# -----------------------------------------------------------------------------
# extracting segments from segs tabledata.
segs = ud.read_upskl(filepath+'segs/'+sname)
# Extracting start and end of segments. Calculating middpoints within periodic box.
U = segs[['U0','U1','U2']].values
V = segs[['V0','V1','V2']].values

#d_skel = ud.nearest_neighbour(pos, segs_mid)[0]
