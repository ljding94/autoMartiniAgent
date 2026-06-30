import os
import numpy as np
from matplotlib import pyplot as plt
import mdtraj as md
from bonded_utilities import *

# Construct a CG model of PEO, load in a AA->CG mapped trajectory, and compute
# bond, angle, and torsion distributions

#----------------#
# Build CG model #
#----------------#

# Number of monomers:
n_mono = 20

# Number of end beads to exclude in bonded distributions
end_mono_exclude = 2 

# Dictionary of beads in monomer:
# Beads are defined by numeric type
bead_data = {}
bead_data[1] = {
    'name': 'SC1',
    'mass': 45, # g/mol
    'charge': 0,  # e
}

# Get bead types for all beads in the chain:
bead_types = {}
j = 0
for mono in range(n_mono):
    for b in list(bead_data):
        bead_types[j] = b
        j += 1
  
# Dictionary of bonds in polymer: 
# In fitting procedure we would store bonded potential parameters here   
bond_data = {}
bond_data[1] = {
    'beads': [1,1],
    }

# Define numeric bead type connecting 2 monomers    
monomer_link = {'start':1,'end':1}    
    
# Angle definitions:
angle_data = {}
angle_data[1] = {
    'beads': [1,1,1],
    }
    
# Torsion definitions:
torsion_data = {}
torsion_data[1] = {
    'beads': [1,1,1,1],
    }
  
# Get bonded dicts mapping type to members:  
bond_list, bond_type_list, bond_type_map, bond_type_reverse_map = get_bond_indices(
    bead_data,
    bead_types,
    bond_data,
    n_mono,
    end_mono_exclude=end_mono_exclude,
    monomer_link=monomer_link,
    )  

angle_list, angle_type_list, angle_type_map, angle_type_reverse_map = get_angle_indices(
    angle_data,
    bond_list,
    bead_data,
    bead_types,
    )  
    
torsion_list, torsion_type_list, torsion_type_map, torsion_type_reverse_map = get_torsion_indices(
    torsion_data,
    angle_list,
    bead_data,
    bead_types,
    )  
    
#---------------------------#
# Load in mapped trajectory #
#---------------------------#

xtc_in = '../../derived/PEO20_solu/PEO20_cg.xtc'
gro_in = '../../derived/PEO20_solu/PEO20_cg.gro'

# Load in a list of single frame dcd from NAMD
traj = md.load(
    xtc_in,
    top=gro_in,
    )

n_frames_traj = traj.n_frames
print(f'Analyzing {n_frames_traj} frames')

#----------------------------#
# Compute bond distributions #
#----------------------------#

# Set histogram parameters:
binwidth_nm = 0.001

hist_data_bonds = compute_bond_distributions(traj,bond_type_map,binwidth_nm)
fit_bonds, rmse_bonds, params_bonds = fit_gaussians(hist_data_bonds)
for t,p in params_bonds.items():
    print(f'  bond type {t}:  μ={p[0]:.4f} nm  σ={p[1]:.4f} nm  RMSE={rmse_bonds[t]:.4f}')

# Set plotting parameters:
plotfile_out_bonds = 'bond_hists.pdf'
bonded_type = 'bond'
xlims = {
    1: (0.25,0.40),
    }

plot_bonded_distributions(
    hist_data_bonds,
    plotfile_out_bonds,
    bonded_type=bonded_type,
    xlims=xlims,
    fit_data=fit_bonds,
    rmse_data=rmse_bonds,
    params=params_bonds,
    )
    
#-----------------------------#
# Compute angle distributions #
#-----------------------------#

# Set histogram parameters:
binwidth_deg = 1

hist_data_angles = compute_angle_distributions(traj,angle_type_map,binwidth_deg)
fit_angles, rmse_angles, params_angles = fit_gaussians(hist_data_angles)
for t,p in params_angles.items():
    print(f'  angle type {t}: μ={p[0]:.2f}°  σ={p[1]:.2f}°  RMSE={rmse_angles[t]:.4f}')

# Set plotting parameters:
plotfile_out_angles = 'angle_hists.pdf'
bonded_type = 'angle'
xlims = {
    1: (0,180),
    }

plot_bonded_distributions(
    hist_data_angles,
    plotfile_out_angles,
    bonded_type=bonded_type,
    xlims=xlims,
    fit_data=fit_angles,
    rmse_data=rmse_angles,
    params=params_angles,
    )

#-------------------------------#
# Compute torsion distributions #
#-------------------------------#

# Set histogram parameters:
binwidth_deg = 1

hist_data_torsions = compute_torsion_distributions(traj,torsion_type_map,binwidth_deg)
fit_torsions, rmse_torsions, params_torsions = fit_gaussians(hist_data_torsions)
for t,p in params_torsions.items():
    print(f'  torsion type {t}: μ={p[0]:.2f}°  σ={p[1]:.2f}°  RMSE={rmse_torsions[t]:.4f}')

# Set plotting parameters:
plotfile_out_torsions = 'torsion_hists.pdf'
bonded_type = 'torsion'
xlims = {
    1: (180,-180),
    }

plot_bonded_distributions(
    hist_data_torsions,
    plotfile_out_torsions,
    bonded_type=bonded_type,
    xlims=xlims,
    fit_data=fit_torsions,
    rmse_data=rmse_torsions,
    params=params_torsions,
    )
    