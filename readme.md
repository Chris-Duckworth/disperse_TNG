# IllustrisTNG DisPerSE catalogues

Cosmic web catalogues for IllustrisTNG 100 and 300. Returns distances to different cosmic web elements for a queried set of points.

Directory contains entire process to make catalogues (including intermediate files that are used internally for disperse). 

**If you are only interested in output (i.e. distances to subhalos or your own chosen set of positions) go straight to /disperse_TNG/...** I've included a handful of Python scripts in ./lib & ./example_scripts to get you started.

Any questions / notice anything weird?   
cduckastro@gmail.com / v0.1 Feb 4 2020.

## Requirements : 
Python with basic packages : numpy, pandas, scipy

## Output :

## Various tracers (used to reconstruct cosmic web) : 

stel_subhalo - cosmic web reconstructed using all subhalos that contain a minimum **stellar** subhalo mass of $10^{8.5} M_{\odot}$. designed to be roughly comparable to what you might be able to recover from observations. only supplied for relatively high redshift.

dm_subalo - cosmic web reconstructed using all subhalos that contain a minimum **total** subhalo mass of $10^{8} M_{\odot}$. designed to more or less map out total DM distribution for underlying filaments. defined for all redshifts.

%gas_density - cosmic web reconstructed using a downsampling of all gas cells to create a gas density field. defined for all redshifts.

## FAQs :

**DisPerSE!?** 

more info can be found here & linked papers for in-depth info: http://www2.iap.fr/users/sousbie/web/html/indexd41d.html 

**Why have you repeated the periodic cube outside each face of the cube to make these catalogues?**

DisPerSE has an option to handle periodic boundaries (i.e. link filaments etc. that cross periodic boundaries). However it struggles with this when using poorly sample datasets (e.g. observational-like reconstruction using galaxy positions). Instead, we have repeated the box in each direction to ensure the filaments here are connected. Just don't go querying points that fall outside of the periodic cube range and everything should be fine.
