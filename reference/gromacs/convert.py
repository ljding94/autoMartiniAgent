import MDAnalysis as mda

# topology/reference structure
u = mda.Universe("step4.0_minimization.gro", "step5_200.dcd")

with mda.Writer("step5_200.xtc", n_atoms=u.atoms.n_atoms) as W:
    for ts in u.trajectory:
        W.write(u.atoms)

