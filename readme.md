# IllustrisTNG DisPerSE cosmic web catalogues

Cosmic web (CW) catalogues for IllustrisTNG 100 and 300. 
Returns distance to filaments and critical points for groups (friend of friend halos) and subhalos for selected snapshots (i.e. all main snapshots between z=2 to z=0).
Directory contains scripts used to create catalogues (including intermediate steps used internally for disperse).

Any questions / notice anything weird / want other snapshots?

chrisjamesduckworth@gmail.com / v0.2 Sept 17 2020.

---

## Only interested in CW distances to FoFs and subhalos?

They are stored in the hdf5 format here: 
They are about 300 Mb (1 Gb) in size per snapshot for TNG 100 (300), and contain the following columns: 

- group_ID (subhalo_ID): index corresponding to row in standard TNG group catalogues.
- d_minima: distance to nearest minimum critical point (void)
- d_saddle_1: distance to nearest 1-saddle point (i.e. critical point where one dimension is collapsing)
- d_saddle_2: distance to nearest 2-saddle point (i.e. critical point where two dimensions are collapsing)
- d_node: distance to nearest maximum critical point (node)
- d_skel: distance to nearest filament segment (computed here as the distance to the nearest segment midpoint).

--- 

## Interested in working directly with the disperse files or want to query your own CW distances for different objects?

In ./lib and ./example_scripts there are a couple of scripts to quickly...
- Make basic visualisations of the identified cosmic web (./example_scripts/plot_web.ipynb)
- Read in disperse files and find own CW distances for a set of points (./example_scripts/find_upskel_distances.py)

### Requirements : 
Python with basic packages : numpy, pandas, scipy, h5py

--- 

## How was the cosmic web identified?

DisPerSE (Discrete Persistent Structure Extractor) is a geometric ridge extractor which identifies critical points in a density field and the unique integral lines between them (i.e. filaments). More info can be found here (plus linked papers): http://www2.iap.fr/users/sousbie/web/html/indexd41d.html 

In a sentence, disperse uses a set of discrete points (e.g. galaxies) to estimate the density field, and hence, identify the morphological features of the cosmic web. Here we reconstruct the cosmic web using two different sets of tracers:

- stel_subhalo [TNG100 and TNG300] - cosmic web reconstructed using all subhalos that contain a minimum **stellar** subhalo mass of 10^{8.5} Msol. Designed to be roughly comparable to what you might be able to recover from observations. Snapshots: 33, 40, 50, 67, 78, 84, 91, 99. (Other snapshots later than 30 available on demand).

- DM_subalo [TNG100 only currently] - cosmic web reconstructed using all subhalos that contain a minimum **total** subhalo mass of 10^{8.5} Msol (i.e. including DM and gas masses). Designed to recover an unbiased cosmic web, however, is less precise than computing from raw particle distrbutions. Snapshots: 33, 40, 50, 67, 78, 84, 91, 99.

---

## Why have you repeated the periodic cube outside each face of the cube to make these catalogues?

DisPerSE has an option to handle periodic boundaries (i.e. link filaments etc. that cross periodic boundaries). However it struggles with this when using poorly sample datasets (e.g. observational-like reconstruction using galaxy positions). Instead, we have repeated the box in each direction to ensure the filaments here are connected. Just don't go querying points that fall outside of the periodic cube range and everything should be fine.