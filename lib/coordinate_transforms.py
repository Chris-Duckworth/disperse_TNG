'''
coordinate_transforms - functions that convert code/physical units and box wrap.
'''

import numpy as np

# defining planck 15 cosmology. astropy is off for some reason.
HubbleParam = 0.6774
OmegaLambda = 0.6911
Omega0 = 0.3089


def H(z):
    term1 = Omega0*(1+z)**3;
    term2 = (1-Omega0-OmegaLambda)*(1+z)**2;
    hubz = 100 * HubbleParam * np.sqrt(term1+term2+OmegaLambda)
    return hubz


def box_wrap(pos_comoving, box_side_length):
    '''This function takes raw comoving (i.e. code units) coordinates (for a given object)
       and wraps them. 
       The first position in the array supplied will be used to centre the object around 
       the origin. Therefore the returned particles will have been shifted and not in 
       origin position relative to the box.     
    '''
    # roughly shifting object to origin.
    dx = pos_comoving - pos_comoving[0]
    
    # if any part of the object falls outside of the box length (centered around 0), then
    # applying periodicity.
    dx[dx > box_side_length/2] = dx[dx > box_side_length/2] - box_side_length
    dx[dx < -box_side_length/2] = dx[dx < -box_side_length/2] + box_side_length
    return dx + pos_comoving[0]


def code_to_physical(pos_comoving, vel_comoving, z):
    ''' This function accepts the code units (ckpc/h for pos) and transforms them to
        physical units including Hubble flow for vel.'''
    
    pos_physical = pos_comoving * 1 / (1 + z) * 1 / HubbleParam
    vel_peculiar = vel_comoving * 1 / np.sqrt(1 + z)
    vel_physical_total = vel_peculiar + H(z) * pos_physical / 1000
    return pos_physical, vel_physical_total


def physical_to_code(pos_physical, vel_physical_total, z):
    ''' This function accepts physical units for pos and vel and transforms 
        them to comoving coordinates. 
        IMPORTANT: This assumes input is the total velocity (i.e. hubble flow included), 
        this will return in code units (i.e. comoving). 
        Make sure that physical pos input are in kpc.'''

    vel_physical_peculiar = vel_physical_total - H(z) * pos_physical / 1000 
    vel_comoving = vel_physical_peculiar * np.sqrt(1 + z)
    pos_comoving = pos_physical * (1 + z) * HubbleParam
    return pos_comoving, vel_comoving


def transform_relative_to_centre(pos, vel, masses=None, potential=None):
    '''Given a set of particle pos and vel, this function transforms these to be relative
       to the overall distribution. 
       By default this finds the centre of mass pos and vel for the set and then 
       re-defines the coordinate system relative to this.
       i.e. pos relative to CoM and vel have overall CoM motion removed.
       
       Alternatively you can define the centre position to be relative to the particle 
       with minimum potential value. Simply provide the potential values for each particle.
       Velocities are always returned relative to the CoM motion.
              
       Masses parameter is optional and if not provided the function will assume that the
       particles should be weighted equally in CoM calculation.
       
       IMPORTANT: Make sure box-wrapping is applied before doing this for border objects.
       
       Parameters
       ----------
       pos :  ndarray
           (n1, Ndim) cartesian (physical and box wrapped) coordinates. (kpc)
       vel :  ndarray
           (n1, Ndim) cartesian (physical and box wrapped) velocities. (km/s * kpc)
       masses : ndarray (optional)
           (n1) masses for each particle.
       potential : ndarray (optional)
           (n1) potential measures for each particle. if supplied will return positions 
           relative to particle with minimum potential.
           
       Returns
       -------
       rel_pos : numpy.ndarray
           (n1, Ndim) positions defined relative to object centre (either CoM or potential)
       rel_vel : numpy.ndarray
           (n1, Ndim) velocities defined relative to CoM motion of object.
       '''
    
    # Set equal weighting if masses not supplied.
    if masses is None:
        masses = np.ones(pos.shape[0])
        
    # Selecting centre definition.
    if potential is None:
        pos_cen = np.sum(pos * masses[:, np.newaxis], axis=0) / np.sum(masses)
    else:
        pos_cen = pos[np.argmin(potential)]

    # Finding CoM motion for all particles.
    CoM_vel = np.sum(vel * masses[:, np.newaxis], axis=0) / np.sum(masses)     

    return pos - pos_cen, vel - CoM_vel
