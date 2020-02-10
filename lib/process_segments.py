'''
process_segments.py

Given a defined skeleton (critical points [.crits] and segments [.segs]) this script deals 
with periodic boundary conditions and re-connects overlapping segments in periodic space.

Chris J Duckworth cduckastro@gmail.com
'''

import numpy as np
import upskl_dist as ud
import coordinate_transforms as ct
import scipy.spatial


def find_segment_inbox(U, V, box_side_length) :
    '''
    Given a set of segments, this selects those that are defined within the periodic box (i.e. both U and V inside 
    the bounds of box_side_length).
    
    ---
    Input :
    
    U : np.array(3, nseg)
        Complete set of 3D positions of starting points of filament segments.
    
    V : np.array
        Complete set of 3D positions of end points of filament segments.
    
    box_side_length : float/int
        Side length of periodic cube/box in same units as U and V.
    
    ---
    Output :
    
    mask : np.array(boolean)
        boolean mask of length nseg
    '''
    # finding all segments that have any dimension outside of the box bounds & setting to False.
    Uoutbox = ~np.any((U > box_side_length) | (U < 0), axis=1)
    Voutbox = ~np.any((V > box_side_length) | (V < 0), axis=1)
    
    # selecting only those which have both U and V entirely inside box.
    return U[Uoutbox & Voutbox], V[Uoutbox & Voutbox]


def move_segment_inbox(U, V, box_side_length) :
    '''
    Given a set of segments, this ensures that both U and V are defined within the periodic box (i.e. inside the bounds
    of box_side_length).
    
    ---
    Input :
    
    U : np.array
        Complete set of 3D positions of starting points of filament segments.
    
    V : np.array
        Complete set of 3D positions of end points of filament segments.
    
    box_side_length : float/int
        Side length of periodic cube/box in same units as U and V.
    
    ---
    Output :
    
    U : np.array
        Complete set of 3D positions of starting points, now defined 0-box_side_length. 
        
    V : np.array
        Complete set of 3D positions of end points, now defined 0-box_side_length. 
    '''
    # creating copies of U and V.
    Utrans = U
    Vtrans = V
    
    # shifting each U and V points inside the box.
    Utrans[Utrans < 0] += box_side_length
    Utrans[Utrans > box_side_length] -= box_side_length
    
    Vtrans[Vtrans < 0] += box_side_length
    Vtrans[Vtrans > box_side_length] -= box_side_length
    return Utrans, Vtrans


def find_segment_matches(Useg, Vseg, Utree, Vtree, tolerance=1):
    '''
    For a supplied segment with 3D start (Useg) and end (Vseg) positions, this searches in the KDtree of all segment 
    points. This assumes that the segment is in the KDtree and hence ignores closest matches for U-U and V-V.
    In-case of identical segment with opposite direction, this also searches for inverse match (i.e) U-V and V-U.
    
    ---
    Input :
    
    Useg : np.array(3)
        3D position of start point of segment.
        
    Vseg : np.array(3)
        3D position of end point of segment.
        
    Utree : scipy.spatial.cKDtree object
        KDtree made up of all segment start points.
        
    Vtree : scipy.spatial.cKDtree object
        KDtree made up of all segment end points.
        
    ---
    Output :
    
    index : int
        Index referring to matched segment. If no match found within the tolerance then this returns -1.
    '''
    
    # finding all combinations of nearest matches between U and V of chosen segment and all of those in the KDtree.
    UUdist, UUind = np.array(Utree.query(Useg, k=2))[:,1]
    VVdist, VVind = np.array(Vtree.query(Vseg, k=2))[:,1]
    UVdist, UVind = np.array(Utree.query(Vseg, k=1))
    VUdist, VUind = np.array(Vtree.query(Useg, k=1))
    
    # checking if segment has a close match in U-U and V-V.
    if (UUdist + VVdist <= tolerance) & (UUind == VVind):
        #print('Match found for U-U and V-V')
        return int(UUind)
    
    elif (UVdist + VUdist <= tolerance) & (UVind == VUind):
        #print('Match found for U-V and V-U')
        return int(UVind)
    
    else :
        #print('No match found')
        return -1


def remove_repeated_segments(U, V, box_side_length, keep_unmatched=False, tolerance=1):
    '''
    Given a set of segments which have been applied to an extended box (i.e. repeated outside of periodic box bounds),
    its very likely that you have multiple filaments in the same positions (once accounting for periodicity). 
    This script finds those that are essentially overlapping and combines them. 
    
    IMPORTANT : This is only applied to filament segments which pass through the walls of the periodic box.
    
    ---
    Input :
    
    U : np.array
    Complete set of 3D positions of starting points of filament segments.
    
    V : np.array
    Complete set of 3D positions of end points of filament segments.
    
    --- 
    Output : 
    
    Uclip : np.array
    Complete set of 3D positions of starting points of filament segments with repeated (in periodic space) removed.
    
    Vclip : np.array
    Complete set of 3D positions of end points of filament segments with repeated (in periodic space) removed.
    '''
    
    # Creating masks for both U and V to find segments that have any dimension outside of the box.
    # assuming box origin is (0,0,0)
    Ubox = np.any((U >= box_side_length) | (U <= 0), axis=1)
    Vbox = np.any((V >= box_side_length) | (V <= 0), axis=1)
    
    # finding all segments which clip through side of the box. (i.e. either U or V falls outside)
    clippers = (~Ubox & Vbox) | (Ubox & ~Vbox)
    
    # slicing clippers into new array.
    Uclip = U[clippers]
    Vclip = V[clippers]
    
    # selecting all clippers that go through upper bound of box and translating to lower box bound.
    # (to identify repeated segments)
    
    # firstly for those with starting points outside upper bound of the box.
    Vclip[Uclip >= box_side_length] = Vclip[Uclip >= box_side_length] - box_side_length
    Uclip[Uclip >= box_side_length] = Uclip[Uclip >= box_side_length] - box_side_length
    
    # now for those with end points outside the upper bound of the box. 
    # by definition the clippers will only have one point outside of these bounds and hence can't be translated twice.
    Uclip[Vclip >= box_side_length] = Uclip[Vclip >= box_side_length] - box_side_length
    Vclip[Vclip >= box_side_length] = Vclip[Vclip >= box_side_length] - box_side_length
    
    # ------------------------------------------------------------------------------------------------
    # now building cKDtree for all segment start and end points. 
    # make sure you have translated your segments to be in same space.
    Utree = scipy.spatial.cKDTree(Uclip)
    Vtree = scipy.spatial.cKDTree(Vclip)
    
    # running through all segments which clip through bounding box. returning those with matches in periodic space.
    matches = []
    for Useg, Vseg in zip(Uclip, Vclip):
        matches.append(find_segment_matches(Useg, Vseg, Utree, Vtree))
    
    matches = np.array(matches) # I know this isn't pythonic but its almost midnight.
    
    # creating 2D array where each row shows the index and its corresponding match.
    match2d = np.array([np.arange(matches.shape[0]), matches]).T
    
    # sorting each row in the matches and then finding unique rows (pairs) only. 
    unique_matches = np.unique( np.sort(match2d, axis=1), axis=0)
    
    # printing total number of matches and unmatched.
    print(str(matches.shape[0])+' total segments crossing through periodic box boundary. \n'+
          str(unique_matches[unique_matches[:,0] != -1].shape[0])+' unique matches (pairs of segments). '+
          str(2*unique_matches[unique_matches[:,0] != -1].shape[0])+' segments combined. \n'+
          str(unique_matches[unique_matches[:,0] == -1].shape[0])+' unmatched segments.')
    
    # checking conditional to keep/remove unmatched segments.
    if keep_unmatched == True :
        # keeping *all* unique segments that clip through bounding box (including those without matches)
        indices = unique_matches[:,0]
        # translating so that all clipper U and Vs are defined with the bounding box.
        Uoutput, Voutput = move_segment_inbox(Uclip[indices], Vclip[indices], box_side_length)
        
    elif keep_unmatched == False :
        indices = unique_matches[unique_matches[:,0] != -1][:,0]
        # translating so that all clipper U and Vs are defined with the bounding box.
        Uoutput, Voutput = move_segment_inbox(Uclip[indices], Vclip[indices], box_side_length)
        
    else : 
        assert False, "keep_umatched must be boolean!"
        
    # finding all segments that are completely defined within the bounding box and appending to clippers.
    Uinbox, Vinbox = find_segment_inbox(U, V, box_side_length)
    print(str(Uinbox.shape[0])+' segments inside the periodic box.')
    
    # combining clippers and entirely in-box segments and returning.
    return np.concatenate((Uoutput, Uinbox), axis=0), np.concatenate((Voutput, Vinbox), axis=0)
