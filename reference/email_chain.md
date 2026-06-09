Did Claude help in building the repository?
From: Ding, Lijie <dingl1@ornl.gov>
Sent: Friday, May 1, 2026 1:24 AM
To: Carrillo, Jan Michael <carrillojy@ornl.gov>; Walker, Christopher <walkercc@ornl.gov>; Kim, Seonghan <kimsn@ornl.gov>
Subject: Re: SMILES → MARTINI 3 CG Pipeline — Setup, Reproduction & Limitations

Hi Guys,
Thank you for the reference, I started building this under this repo.
https://github.com/ljding94/autoMartiniAgent

Regards,
Lijie

From: Carrillo, Jan Michael <carrillojy@ornl.gov>
Date: Thursday, April 30, 2026 at 9:02 PM
To: Walker, Christopher <walkercc@ornl.gov>; Kim, Seonghan <kimsn@ornl.gov>; Ding, Lijie <dingl1@ornl.gov>
Subject: Re: SMILES → MARTINI 3 CG Pipeline — Setup, Reproduction & Limitations

Can you guys take a look at this as well
https://doi.org/10.1021/acs.jcim.5c02903
From: Walker, Christopher <walkercc@ornl.gov>
Sent: Wednesday, April 29, 2026 1:14:45 PM
To: Kim, Seonghan <kimsn@ornl.gov>; Ding, Lijie <dingl1@ornl.gov>
Cc: Carrillo, Jan Michael <carrillojy@ornl.gov>
Subject: RE: SMILES → MARTINI 3 CG Pipeline — Setup, Reproduction & Limitations

Thanks Seonghan,

I think this workflow looks good.

I tried to use the PSBMA monomer SMILES string (generated from LigParGen) with automartini:

python -m auto_martini --smi "CC(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]" --mol psbma --cg psbma_monomer_cg.gro --top psbma_monomer.itp -–fpred

It did not converge in over 30 minutes and seems stuck.

The –-fpred flag is needed, otherwise it will fail:
ALOGPS can't predict fragment: [H]C([H])([H])C([H])(C(=O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]

A smaller polymer sidechain PMETAC (“CC(C)C(=O)OCC[N+](C)(C)”) does succeed in ~1 minute, but the mapping doesn’t follow basic Martini principles, putting 2-3 heavy atoms in standard size beads which should get 4  – it gives:

[atoms]
; id    type    resnr   residue  atom    cgnr    charge  smiles
    1     C3      1     pmetac     C01     1         0   ; CCC
    2     P1      1     pmetac     P01     2         0   ; OC=O
    3     C5      1     pmetac     C02     3         0   ; CC
    4     Qd      1     pmetac     Q01     4         1   ; C[N+]C

This should be fine for a starting point, but I think we want the agent to be able to manipulate the number of CG beads in the mapping and change Martini sizes as needed.

Chris

From: Kim, Seonghan <kimsn@ornl.gov>
Sent: Wednesday, April 29, 2026 11:21 AM
To: Ding, Lijie <dingl1@ornl.gov>
Cc: Carrillo, Jan Michael <carrillojy@ornl.gov>; Walker, Christopher <walkercc@ornl.gov>
Subject: FW: SMILES → MARTINI 3 CG Pipeline — Setup, Reproduction & Limitations

Hi All,

I have attached a workflow for handling small molecules using AutoMARTINI3 below. The details can be certainly adjusted but I think this could be a good starting point. Please let me know if you have any questions or suggestions.

For reference, I also included the detailed AutoMARTINI3/GROMACS reproduction notes below that I shared with Lijie.

Thanks,
Seonghan

----

Stage 1. Initial CG generation
Input SMILES / molecule type
    ↓
Agent classifier:
small molecule ≤ 25 heavy atoms
    → AutoMARTINI3
medium molecule / oligomer
    → fragment + AutoMARTINI3 + assembly agent
polymer / repeating unit
    → Polyply agent
protein / peptide
    → Martinize2 agent

Stage 2: Refinement
atomistic reference trajectory
    → bonded distribution target generation
    → BI or Swarm-CG
    → refined .itp

Stage 3: Validation
CG simulation
    → Comparison with AA-mapped distribution
    → report


From: Kim, Seonghan <kimsn@ornl.gov>
Date: Wednesday, April 29, 2026 at 12:55 AM
To: Ding, Lijie <dingl1@ornl.gov>
Subject: SMILES → MARTINI 3 CG Pipeline — Setup, Reproduction & Limitations

Hi Lijie,

Here's a guide to reproduce the SMILES → MARTINI 3 coarse-grained MD pipeline. The goal is to go from a SMILES string to an energy-minimized CG system using Auto-MartiniM3 and GROMACS, which I think AI Agent should be able to handle.

---

REQUIREMENTS

- GROMACS >= 2021
- conda
- git

---

SETUP

1. Create conda environment:
   $conda create -n autom3 python=3.10 -y
 $conda activate autom3

2. Install Auto-MartiniM3:
 $git clone https://github.com/Martini-Force-Field-Initiative/Automartini_M3.git
 $cd Automartini_M3
 $pip install -e .
 $pip install rdkit
 $cd ..

3. Download MARTINI 3 force field files into ff/ directory:
 $mkdir -p ff && cd ff
 $wget https://cgmartini-library.s3.ca-central- 1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_v3.0.0.itp
 $wget https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_v3.0.0_solvents_v1.itp
 $wget https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/example_applications/solvent_systems/water.gro
 $cd ..

---

DIRECTORY STRUCTURE

   cg_pipeline/
   ├── ff/
   │   ├── martini_v3.0.0.itp
   │   ├── martini_v3.0.0_solvents_v1.itp
   │   └── water.gro
   └── test_octanol/
       ├── em.mdp
       ├── topol.top
       └── (run pipeline here)

---

NOTE: MOL NAME CONFLICT

Do NOT use MOL="OCT" — it conflicts with octane already defined in martini_v3.0.0_solvents_v1.itp. Use "OCOL" for 1-octanol instead. In general, check for name conflicts against the solvents itp before running.

---

RUN SCRIPT (save as run.sh inside test_octanol/)

   #!/usr/bin/env bash
   set -euo pipefail

   SMI="CCCCCCCCO"
   MOL="OCOL"
   NMOL=50
   BOX=4

   # 1. Generate CG topology from SMILES
   $python -m auto_martiniM3 --smi "$SMI" --mol "$MOL" --canon -v

   # 2. Pack molecules into box
   $gmx insert-molecules -ci ${MOL}.gro -nmol $NMOL -box $BOX $BOX $BOX -o packed.gro

   # 3. Solvate with MARTINI 3 water
   $gmx solvate -cp packed.gro -cs ../ff/water.gro -o solv.gro -p topol.top -radius 0.21

   # 4. Prepare EM input
   $gmx grompp -f em.mdp -c solv.gro -p topol.top -o em.tpr

   # 5. Run energy minimization
   $gmx mdrun -deffnm em -v

   $echo "=== EM converged ==="
   $tail -5 em.log

---

topol.top (create manually before running):

   #include "../ff/martini_v3.0.0.itp"
   #include "../ff/martini_v3.0.0_solvents_v1.itp"
   #include "OCOL.itp"

   [ system ]
   Octanol in MARTINI 3 water

   [ molecules ]
   OCOL 50

Note: gmx solvate will automatically append the water (W) count to [ molecules ].

---

em.mdp (MARTINI 3 standard settings; energy minimization):

   integrator        = steep
   nsteps            = 5000
   emtol             = 100
   emstep            = 0.01
   cutoff-scheme     = Verlet
   nstlist           = 20
   coulombtype       = reaction-field
   rcoulomb          = 1.1
   epsilon_r         = 15
   vdw-type          = cutoff
   vdw-modifier      = Potential-shift-verlet
   rvdw              = 1.1
   pbc               = xyz
   constraints       = none

---

WHAT TO CHECK AFTER RUNNING

1. OCOL.itp exists and contains [ bonds ] and [ angles ] sections (not just [ atoms ])
2. gmx grompp passes with no ERRORs
3. em.gro is generated after mdrun
4. Potential Energy in em.log is finite and negative (~-10^4 to -10^5 kJ/mol)
5. Maximum force converges below emtol (100 kJ/mol/nm)

Quick check:
   grep 'Potential Energy\|Maximum force' em.log | tail -5

---

LIMITATIONS OF THIS APPROACH

Auto-MartiniM3 is suitable for small organic molecules but has the following constraints:

- Hard limit of 25 heavy atoms per molecule. Larger drug-like molecules are not supported.
- For polymers, Polyply (https://github.com/marrink-lab/polyply_1.0) is the appropriate tool — it takes a monomer .itp and generates arbitrary-length polymer topologies and coordinates. Stock MARTINI 3 library covers PEO, PS, and a few others.
- No built-in validation. The pipeline does not verify whether generated parameters reproduce experimental observables (density, partition free energy, etc.).

---

Let me know if you run into any issues or have things to discuss.

Best,
Seonghan
