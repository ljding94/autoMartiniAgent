import os
import numpy as np
import mdtraj as md
from matplotlib import pyplot as plt

# Tools for getting bond, angle, and torsion lists from a CG model,
# and computing, plotting bonded histograms.

def get_bond_indices(bead_data,bead_types,bond_data,n_mono,end_mono_exclude=0,monomer_link={'start':1,'end':1}):
    """
    Given bead definitions (bead_data) and bond_data defining which bead types are
    connected by bonds, determine the bead indices grouped by their
    bond type. This is used for computing bond distributions.
    
    ***Note - this function is for polymers with linear CG side chains only (linear or PMMA-like)
    
    :param bead_data: dictionary mapping numeric bead type index to name, mass, charge
    :type bead_data: dict(1: dict('name':str,'mass':float,'charge':float))
    
    :param bond_data: dictionary mapping numeric bond type to 1) numeric atom types involved in the bond, 2) force constant (kJ/mol/nm^2), 3) eq. distance (nm), 4) gromacs function type
    :type bond_data: dict(1: dict('beads':list(int),'kb':float,'b0':float,'gromacs_func':int)    
    
    :param end_mono_exclude: exclude particles from this many monomers from each chain end (default=0)
    :type end_mono_exclude: int    
   
    :param monomer_link: dictionary containing the bead types which are bonded between monomers (default={'start':1,'end':1})
    :type monomer_link: dict('start':int,'end':int)
   
    returns:
        bond_list - list of bead index pairs bonded together in the chain ( list([int,int]) )
        bond_type_list - bond types corresponding the entries in bond_list ( list([int,int]) )
        bond_type_map - dictionary mapping bond numeric type to bond pairs of that type ( dict{int:(list([int,int]))} )
        bond_type_reverse_map - reverse lookup dictionary mapping pair of bead numeric types to the bond numeric type ( dict{[int,int]:int} )     
    """    
    
    n_bead_per_mono = len(bead_data)
    n_bead_per_chain = n_mono*n_bead_per_mono

    # Map numeric bead types to bond types: 
    bond_type_reverse_map = {}
    bond_type_map = {}  # Maps bond type to list of all relevant bonds over all chains    
    
    for key,val in bond_data.items():
        bond_type_map[key] = []    
    
    for key,val in bond_data.items():
        # Here key is numeric bond type, val is the inner dictionary
        beads = val['beads']            # For example, [1,2]
        bond_type_reverse_map[str(beads)] = key # For example, map '[1,2]' to bond type 1
        bond_type_reverse_map[str(list(reversed(beads)))] = key # Reverse sequence maps to same

    # Get the monomer-monomer linking bond type:
    linker_bead0_type = monomer_link['start']
    linker_bead1_type = monomer_link['end']
    linker_bond = [linker_bead0_type,linker_bead1_type]
    
    linker_bond_type = bond_type_reverse_map[str(linker_bond)]

    bond_list = []
    bond_type_list = []

    for mono in range(end_mono_exclude,n_mono-end_mono_exclude):
        # Add intra-monomer bonds (anything not of type monomer_link)
        
        # Get the bead indices involved in current monomer:
        mono_beads = {} # actual indices of the beads
        template_beads = {} # indices mapped to the template
        for b in range(n_bead_per_mono):
            mono_beads[b+1] = mono*n_bead_per_mono + b + 1
            
        monomer_bond_list = []
        for pair_str in bond_type_reverse_map:
            # Note - bond_type_reverse_map also contains the reversed sequence,
            # so need to check if the bond has already been added.
            
            # These bead pairs are stored as strings - convert back to a list:
            pair = list(eval(pair_str[1:-1]))
            
            # Bead indices:
            bead0 = mono_beads[pair[0]]
            bead1 = mono_beads[pair[1]]
            
            # Bead types:
            bead0_type = bead_types[bead0]
            bead1_type = bead_types[bead1]
            
            # Here we check if the forward and reverse is already in the bond list,
            # and check if the bond involves the same bead type (such as SC1-SC1 backbone links)
            if [bead0,bead1] not in monomer_bond_list and [bead1,bead0] not in monomer_bond_list and bead0 != bead1:
                monomer_bond_list.append([bead0,bead1])

        for mono_bond in monomer_bond_list:
            bead_types_curr = [
                bead_types[mono_bond[0]],
                bead_types[mono_bond[1]],
                ]
            bond_type_curr = bond_type_reverse_map[str(bead_types_curr)]
            
            # For polystyrene, the 1-2 bond is both an intra-monomer type and a linker type
            # For linear PDMAEMA type Martini models, the 1-1 bond is the linker type and shouldn't be added here
            
            if linker_bond_type != [1,1]:
                bond_list.append(mono_bond)
                bond_type_list.append(bond_type_curr)
                bond_type_map[bond_type_curr].append(mono_bond)
            
        for mono in range(end_mono_exclude,n_mono-1-end_mono_exclude):
            # Add the monomer-monomer linker bonds:
            bead1 = mono*n_bead_per_mono + monomer_link['start']
            bead1next = (mono+1)*n_bead_per_mono + monomer_link['end']
            
            bond_indices_curr = [bead1,bead1next]
            bead_types_curr = [
                bead_types[bead1],
                bead_types[bead1next]
                ]
            bond_type_curr = bond_type_reverse_map[str(bead_types_curr)]
            
            bond_list.append(bond_indices_curr)
            bond_type_list.append(bond_type_curr) 
            bond_type_map[bond_type_curr].append(bond_indices_curr)         
        
    return bond_list, bond_type_list, bond_type_map, bond_type_reverse_map
    
    
def get_angle_indices(angle_data,bond_list,bead_data,bead_types):
    """
    Given bead definitions (bead_data) and angle_data defining which bead types make
    up each angle, determine the bead indices grouped by their
    angle type. This is used for computing angle distributions.
    
    :param angle_data: dictionary mapping numeric angle type to 1) numeric atom types involved in the angle, 2) force constant (kJ/mol/nm^2), 3) eq. distance (nm), 4) gromacs function type
    :type angle_data: dict(1: dict('beads':list(int),'ka':float,'theta0':float,'gromacs_func':int)
    
    :param bond_list: list bead index pairs bonded together
    :type bond_list: list(list([int,int]))
    
    :param bead_data: dictionary mapping numeric bead type index to name, mass, charge
    :type bead_data: dict(1: dict('name':str,'mass':float,'charge':float))
    
    returns:
        angle_list - list of bead index triplets making up the angles ( dict{int:(list([int,int,int]))} )
        angle_type_list - angle types corresponding the entries in angle_list ( dict{int:(list([int,int,int]))} )
        angle_type_map - dictionary mapping angle numeric type to angle triplets of that type across all chains ( dict{int:(list([int,int,int]))} )
        angle_type_reverse_map - reverse lookup dictionary mapping triplet of bead numeric types to the angle numeric type ( dict{[int,int,int]:int} )     
    """    
    
    angle_list = {}
    angle_type_list = {}
    angle_type_map = {}  # Maps angle type to list of all relevant angles
    
    for key,val in angle_data.items():
        angle_type_map[key] = []    
       
    # Map numeric bead types to angle types: 
    angle_type_reverse_map = {}
    
    for key,val in angle_data.items():
        # Here key is numeric angle type, val is the inner dictionary
        beads = val['beads']            # For example, [1,2,3]
        angle_type_reverse_map[str(beads)] = key # For example, map '[1,2,3]' to angle type 1
        angle_type_reverse_map[str(list(reversed(beads)))] = key # Reverse sequence maps to same

    angle_list = []
    angle_type_list = []
    
    for bond1 in bond_list:
        for bond2 in bond_list:
            # Check if bonds are different:
            if bond1 != bond2 and bond1[::-1] != bond2:
                # Make sure the bonds share a common atom
                bond_angle = []
                if bond2[0] in bond1 or bond2[1] in bond1:
                    if bond2[0] == bond1[1]:
                        bond_angle = [bond1[0], bond1[1], bond2[1]]
                    elif bond2[0] == bond1[0]:
                        bond_angle = [bond1[1], bond1[0], bond2[1]]
                    elif bond2[1] == bond1[1]:
                        bond_angle = [bond1[0], bond1[1], bond2[0]]
                    elif bond2[1] == bond1[0]:
                        bond_angle = [bond1[1], bond1[0], bond2[0]]
                if len(bond_angle) > 0:
                    # Check if already in angle_list and add:
                    if bond_angle not in angle_list and bond_angle[::-1] not in angle_list:
                        angle_list.append(bond_angle)
                        bead_type1 = bead_types[bond_angle[0]]
                        bead_type2 = bead_types[bond_angle[1]]
                        bead_type3 = bead_types[bond_angle[2]]
                        
                        # Get angle type from bead types:
                        angle_bead_types = [bead_type1,bead_type2,bead_type3]
                        angle_type_curr = angle_type_reverse_map[str(angle_bead_types)]
                        angle_type_list.append(angle_type_curr)
                        
                        # Add angle to type-specific dictionary:
                        angle_type_map[angle_type_curr].append(bond_angle)
                    
    return angle_list, angle_type_list, angle_type_map, angle_type_reverse_map        
    

def get_torsion_indices(torsion_data,angle_list,bead_data,bead_types):
    """
    Given bead definitions (bead_data) and torsion_data defining which bead types make
    up each torsion, determine the bead indices grouped by their
    torsion type. This is used for computing torsion distributions.
    
    :param torsion_data: dictionary mapping numeric torsion type to numeric atom types involved in the torsion
    :type torsion_data: dict(int: dict('beads':list(int)))
    
    :param angle_list: list bead index triplets making up each angle
    :type angle_list: list(list([int,int,int]))
    
    :param bead_data: dictionary mapping numeric bead type index to name, mass, charge
    :type bead_data: dict(1: dict('name':str,'mass':float,'charge':float))  
      
    returns:
        torsion_list - list of bead index quartets making up the torsions ( list([int,int,int,int]) )
        torsion_type_list - torsion types corresponding the entries in torsion_list ( list([int,int,int,int]) )
        torsion_type_map - dictionary mapping torsion numeric type to torsion quartets of that type across all chains ( dict{int:(list([int,int,int,int]))} )
        torsion_type_reverse_map - reverse lookup dictionary mapping quartet of bead numeric types to the torsion numeric type ( dict{[int,int,int,int]:int} )     
    """      
    
    torsion_list = {}
    torsion_type_list = {}
    torsion_type_map = {}  # Maps torsion type to list of all relevant torsions over all chains
    
    for key,val in torsion_data.items():
        torsion_type_map[key] = []    
       
    # Map numeric bead types to torsion types: 
    torsion_type_reverse_map = {}
    
    for key,val in torsion_data.items():
        # Here key is numeric torsion type, val is the inner dictionary
        beads = val['beads']               # For example, [1,1,2,3]
        torsion_type_reverse_map[str(beads)] = key # For example, map '[1,1,2,3]' to torsion type 1
        torsion_type_reverse_map[str(list(reversed(beads)))] = key # Reverse sequence maps to same

    # Search for torsions as 2 overlapping angles with a shared bond:
    torsion_list = []
    torsion_type_list = []
    
    for i in range(len(angle_list)):
        for j in range(i+1,len(angle_list)):
            angle_1 = angle_list[i]
            angle_2 = angle_list[j]

            torsion = []

            # Check overlap:
            if [angle_1[1],angle_1[2]] == [angle_2[0],angle_2[1]]:
                # 0 1 2
                #   0 1 2
                torsion = [angle_1[0], angle_1[1], angle_1[2], angle_2[2]]
                    
            elif [angle_1[1],angle_1[2]] == [angle_2[2],angle_2[1]]:
                # 0 1 2
                #   2 1 0
                torsion = [angle_1[0], angle_1[1], angle_1[2], angle_2[0]]
                    
            elif [angle_1[1],angle_1[0]] == [angle_2[0],angle_2[1]]:
                # 2 1 0
                #   0 1 2
                torsion = [angle_1[2], angle_1[1], angle_1[0], angle_2[2]]
                    
            elif [angle_1[1],angle_1[0]] == [angle_2[2],angle_2[1]]:
                # 2 1 0
                #   2 1 0
                torsion = [angle_1[2], angle_1[1], angle_1[0], angle_2[0]]
            
            if len(torsion) > 0:
                # Check that torsion is new and that ends are not the same particle:
                if torsion not in torsion_list and torsion[::-1] not in torsion_list and torsion[0] != torsion[3]:
                    torsion_list.append(torsion)

                    bead_type1 = bead_types[torsion[0]]
                    bead_type2 = bead_types[torsion[1]]
                    bead_type3 = bead_types[torsion[2]]
                    bead_type4 = bead_types[torsion[3]]
                    
                    # Get angle type from bead types:
                    torsion_bead_types = [bead_type1,bead_type2,bead_type3,bead_type4]
                    torsion_type_curr = torsion_type_reverse_map[str(torsion_bead_types)]
                    torsion_type_list.append(torsion_type_curr)
                    
                    # Add angle to type-specific dictionary:
                    torsion_type_map[torsion_type_curr].append(torsion)

    return torsion_list, torsion_type_list, torsion_type_map, torsion_type_reverse_map    


def compute_bond_distributions(traj,bond_type_map,binwidth_nm=0.001):
    """
    Compute bond stretching histograms from mdtraj traj object using bond definitions from bond_type_map
    
    :param traj: mdtraj trajectory
    :type traj: mdtraj.Trajectory
    
    :param bond_type_map: dictionary mapping numeric bond type to list of member pairs
    :type bond_type_map: dict{int:list(int,int)}
    
    :param binwidth_nm: histogram bin width in nm (default=0.001)
    :type binwidth_nm: float
    
    returns:
        hist_data_bonds - dict{bond_type:2d numpy.array}
    """
    bin_edges = np.arange(0.0,0.8+binwidth_nm,binwidth_nm)
    bin_centers = np.zeros((len(bin_edges)-1))
    for i in range(len(bin_edges)-1):
        bin_centers[i] = (bin_edges[i]+bin_edges[i+1])/2

    # Initialize dict for bond histogram data
    hist_data_bonds = {}
    bond_types = list(bond_type_map)
    
    # Compute the bond distributions:
    for b,bond_list in bond_type_map.items():
        # This returns an [nframes x n_bonds] array:
        bonded_val_array = md.compute_distances(
            traj,
            bond_list,
            opt=True,
            )
            
        # Reshape as a 1D array for histogramming:
        array_shape = bonded_val_array.shape  
        bonded_val_array = np.reshape(bonded_val_array,(array_shape[0]*array_shape[1],1))

        # Histogram the data:    
        h_out, bin_edges_out = np.histogram(
            bonded_val_array,
            bins=bin_edges,
            density=True,
            )

        array_out = np.zeros((2,len(h_out)))
        array_out[0,:] = bin_centers
        array_out[1,:] = h_out

        hist_data_bonds[b] = array_out
        
    return hist_data_bonds    
    
    
def compute_angle_distributions(traj,angle_type_map,binwidth_deg):
    """
    Compute bond angle bending histograms from mdtraj traj object using angle definitions from angle_type_map
        
    :param traj: mdtraj trajectory
    :type traj: mdtraj.Trajectory
    
    :param angle_type_map: dictionary mapping numeric angle type to list of member triplets
    :type angle_type_map: dict{int:list(int,int,int)}
    
    :param binwidth_deg: histogram bin width in degrees (default=1)
    :type binwidth_deg: float
        
    returns:
        hist_data_angles - dict{angle_type:2d numpy.array}
    """
    bin_edges = np.arange(0,180+binwidth_deg,binwidth_deg)
    bin_centers = np.zeros((len(bin_edges)-1))
    for i in range(len(bin_edges)-1):
        bin_centers[i] = (bin_edges[i]+bin_edges[i+1])/2

    # Initialize dict for angle histogram data
    hist_data_angles = {}
    angle_types = list(angle_type_map)
    
    # Compute the angle distributions:
    for a,angle_list in angle_type_map.items():
        # This returns an [nframes x n_angles] array:
        bonded_val_array = md.compute_angles(
            traj,
            angle_list,
            opt=True,
            )
            
        # Convert to degrees and reshape as a 1D array for histogramming:
        array_shape = bonded_val_array.shape  
        bonded_val_array = (180/np.pi)*np.reshape(bonded_val_array,(array_shape[0]*array_shape[1],1))
        
        # Histogram the data:    
        h_out, bin_edges_out = np.histogram(
            bonded_val_array,
            bins=bin_edges,
            density=True,
            )

        array_out = np.zeros((2,len(h_out)))
        array_out[0,:] = bin_centers
        array_out[1,:] = h_out

        hist_data_angles[a] = array_out
        
    return hist_data_angles        
    
    
def compute_torsion_distributions(traj,torsion_type_map,binwidth_deg):
    """
    Compute torsion histograms from mdtraj traj object using torsion definitions from torsion_type_map
        
    :param traj: mdtraj trajectory
    :type traj: mdtraj.Trajectory
    
    :param torsion_type_map: dictionary mapping numeric torsion type to list of member quartets
    :type torsion_type_map: dict{int:list(int,int,int,int)}
    
    :param binwidth_deg: histogram bin width in degrees (default=1)
    :type binwidth_deg: float        
        
    returns:
        hist_data_torsions - dict{torsion_type:2d numpy.array}
    """
    bin_edges = np.arange(-180,180+binwidth_deg,binwidth_deg)
    bin_centers = np.zeros((len(bin_edges)-1))
    for i in range(len(bin_edges)-1):
        bin_centers[i] = (bin_edges[i]+bin_edges[i+1])/2

    # Initialize dict for torsion histogram data
    hist_data_torsions = {}
    torsion_types = list(torsion_type_map)
    
    # Compute the angle distributions:
    for t,torsion_list in torsion_type_map.items():
        # This returns an [nframes x n_angles] array:
        bonded_val_array = md.compute_dihedrals(
            traj,
            torsion_list,
            opt=True,
            )
            
        # Convert to degrees and reshape as a 1D array for histogramming:
        array_shape = bonded_val_array.shape  
        bonded_val_array = (180/np.pi)*np.reshape(bonded_val_array,(array_shape[0]*array_shape[1],1))
        
        # Histogram the data:    
        h_out, bin_edges_out = np.histogram(
            bonded_val_array,
            bins=bin_edges,
            density=True,
            )

        array_out = np.zeros((2,len(h_out)))
        array_out[0,:] = bin_centers
        array_out[1,:] = h_out

        hist_data_torsions[t] = array_out
        
    return hist_data_torsions     
    

def plot_bonded_distributions(hist_data,plotfile_out,bonded_type='bond',xlims={},ylims={},ncol=1):
    """
    Plot bonded distributions with each type on a separate subplot
    
    :param data_bonds_multi: dictionary containing bonded histogram data organized by numeric type
    :type density_data: dict{int:2d numpy array}       
    
    :param plotfile_out: file name for output plot 
    :type plotfile_out: str    
    
    :param bonded_type: type of bonded distribution being plotted ('bond','angle','torsion')
    :type bonded_type: str    
    
    :param xlims: map numeric type to x limits
    :type xlims: dict(int:tuple)
    
    :param ylims: map numeric type to y limits
    :type ylims: dict(int:tuple)       
    
    :param ncol: number of columns in the subplot structure (default=2)
    :type ncol: int    
    
    """
    bond_type_list = list(hist_data)
    n_bonds = len(bond_type_list)
    
    # Determine optimal subplot layout with specified number of columns
    nrow = int(np.ceil(n_bonds/ncol))
    
    # Set the calls to axes (can be 0d, 1d, or 2d based on the number of bond types):
    if nrow == 1 and ncol == 1:
        ax_call = 'ax'
    elif nrow == 1 and ncol > 1:
        ax_call = 'ax[col]'
    elif nrow > 1 and ncol == 1:
        ax_call = 'ax[row]'
    elif nrow > 1 and ncol > 1:
        ax_call = 'ax[row,col]'      
        
    # Set figure size:
    if nrow > 2:
        figsize = (10,3*nrow)
    else:
        figsize = (10,6)
        
    fig,ax = plt.subplots(nrows=nrow,ncols=ncol,figsize=figsize)        
    
    subplot_id = 0
        
    for i,bond in enumerate(bond_type_list): 
        # Determine subplot index:
        row = int(np.floor(subplot_id/(ncol)))
        col = int(subplot_id%ncol)    
            
        eval(ax_call).plot(
            hist_data[bond][0,:],
            hist_data[bond][1,:],
            )
            
        eval(ax_call).set_title(f'{bonded_type} type {bond}')
    
        subplot_id += 1
        
    # Format the subplots:
    if bonded_type == 'bond':
        xlabel_str = 'Bond stretching distance (nm)'
    elif bonded_type == 'angle':
        xlabel_str = 'Bending angle (degrees)'
    elif bonded_type == 'torsion':
        xlabel_str = 'Torsion angle (degrees)'
        
    bond_id = 0
    for row in range(nrow):
        for col in range(ncol):
            bond_name = bond_type_list[bond_id]
            
            ylim = eval(ax_call).get_ylim()
            eval(ax_call).set_ylim([0,ylim[1]])
            
            if bond_name in xlims:
                # Set specified x limits:
                eval(ax_call).set_xlim(xlims[bond_name])
            
            if bond_name in ylims:
                # Set specified y limits:
                eval(ax_call).set_ylim(ylims[bond_name])

            bond_id += 1
                
    # Add the x and y labels at edge subplots:            
    for row in range(nrow):
        for col in range(ncol):
            if row == nrow-1:
                eval(ax_call).set_xlabel(xlabel_str,fontsize=14)
            if col == 0:
                eval(ax_call).set_ylabel('Probability',fontsize=14)

    plt.subplots_adjust(
        left=0.1,
        right=0.9,
        hspace=0.3,
        wspace=0.15,
        top=0.9,
        bottom=0.1,
        )

    fig.savefig(plotfile_out)  
    plt.close()        