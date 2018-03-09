#!/usr/bin/env python

# Script to run by-residue SASA calculations for LeuT.
# For speed, the system is truncated to only include the
# closest 100 lipid residues to the protein using pytraj.
# MDtraj is used for the SASA calculation. 
# Relative ASA is normalised to GXG areas from
# Miller S. et. al., J. Mol. Biol., 1987, 196(3), 641-656

import pytraj as pt
import mdtraj as md
import numpy as np
from datetime import datetime
import copy

relatives = { 'A' : 113., 'R' : 241., 'N' : 158., 'D' : 151.,\
              'C' : 140., 'Q' : 189., 'E' : 183., 'G' : 85., \
              'H' : 194., 'I' : 182., 'L' : 180., 'K' : 211., \
              'M' : 204., 'F' : 218., 'P' : 143., 'S' : 122., \
              'T' : 146., 'W' : 259., 'Y' : 229., 'V' : 160., \
              None : 1.} # For lipids without a code

# Both parmfiles must already be available. Won't be created here.
orig_trajn = '../2us_full_run2.dcd' 
orig_parmn = '../../3tt3_y268a_full.psf'

trimmed_trajn = '3TT3_100lipids_tmp.dcd'
trimmed_parmn = '../../3TT3_100_lipids.psf'

#frames = range(0,10) # None = all frames
#Set up outputs
#with open("Protein_SASA_byframe.dat", 'w') as f:
#    f.write("# Residue ASA / A^2. Rows = residues, columns = frames\n")
#with open("Protein_relSASA_byframe.dat", 'w') as f:
#    f.write("# Residue ASA, normalised to G-X-G area. Rows = residues, columns = frames\n")

# First/last residue = first protein, last lipid (1-indexed)
first_protein = 1
last_protein = 507
first_lipid = 508
last_lipid = 755
all_prot_sasas = np.empty((0,507)) # length of protein residues
all_rel_prot_sasas = np.empty((0,507)) # length of protein residues
# Chunk trajectories
for endframe in range(1000,21000,1000):
    startframe = endframe - 1000
    frames = range(startframe,endframe)
    # First/last residue = first protein, last lipid (1-indexed)
    first_protein = 1
    last_protein = 507
    first_lipid = 508
    last_lipid = 755

    # Read & trim to closest 100 with pytraj
    pt_traj = pt.load(orig_trajn, top=orig_parmn, mask=':%d-%d,%d-%d' \
                      % (first_protein, last_protein, first_lipid, last_lipid),\
                      frame_indices=frames)
    print "Read in pytraj to frame %d" % endframe
    print datetime.now()

    trimmed_pt_traj = pt.closest(pt_traj, mask=':%d-%d' % (first_protein, last_protein),\
                                 solvent_mask=':%d-%d' % (first_lipid, last_lipid),\
                                 n_solvents=100, dtype='trajectory')
    print "Done closest to frame %d" % endframe
    print datetime.now()
    pt.write_traj(trimmed_trajn, trimmed_pt_traj, overwrite=True)

    # MDtraj SASA calculation
    md_traj = md.load(trimmed_trajn, top=trimmed_parmn)
    all_sasas = md.shrake_rupley(md_traj, mode='residue')
    print "Done SASA to frame %d" % endframe
    print datetime.now()
    all_sasas *= 100 # Convert to A^2
    prot_sasas = all_sasas[:,first_protein-1:last_protein]
    rel_prot_sasas = copy.deepcopy(prot_sasas)
    for idx, res in enumerate(md_traj.top.subset(md_traj.top.select("resid %d to %d"\
                              % (first_protein-1, last_protein-1))).residues):
        rel_prot_sasas[:,idx] /= relatives[res.code]
    print "Done relSASA to frame %d" % endframe
    print datetime.now()
    all_prot_sasas = np.append(all_prot_sasas, prot_sasas, axis=0)
    all_rel_prot_sasas = np.append(all_rel_prot_sasas, rel_prot_sasas, axis=0)
# save
np.savetxt("Protein_SASA_byframe.dat", all_prot_sasas.T, header="Residue ASA / A^2. Rows = residues, columns = frames", fmt='%8.3f')
np.savetxt("Protein_relSASA_byframe.dat", all_rel_prot_sasas.T, header="Residue ASA, normalised to G-X-G area. Rows = residues, columns = frames", fmt='%8.5f')
# Matrix casting & transpose reshapes into row-wise residues
np.savetxt("Protein_SASA_mean.dat", np.matrix(np.mean(all_prot_sasas, axis=0)).T, header="Residue ASA / A^2", fmt='%8.3f')
np.savetxt("Protein_relSASA_mean.dat", np.matrix(np.mean(all_rel_prot_sasas, axis=0)).T, header="Residue ASA, normalised to G-X-G area", fmt='%8.5f')



