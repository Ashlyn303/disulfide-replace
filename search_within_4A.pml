# Check for input file availability
python
import os
input_pdb = "input_files/5B3N.cif"
if not os.path.exists(input_pdb):
    print(f"ERROR: Input file {input_pdb} not found. Stopping.")
    quit()
python end

# Load the specific structure
load input_files/5B3N.cif

# Select disulfide Cys sidechain atoms (CB, SG)
select ds_atoms, (resn CYS and (name CB or name SG))

# Find L, I, M neighbors with sidechain atoms within 4A of disulfide sidechains
# Excludes backbone N, CA, C, O
select near_lim, (byres (resn LEU+ILE+MET and not (name N+CA+C+O)) within 4 of ds_atoms)

# Print specific residue info for the combinatorial script
python
from pymol import stored
stored.lim_info = []
stored.cys_info = []

# Collect (Chain, Name, Index) for neighbors and disulfides
cmd.iterate("(near_lim and name CA)", "stored.lim_info.append((chain, resn, resi))")
cmd.iterate("(resn CYS and name CA)", "stored.cys_info.append((chain, resn, resi))")

print("\n--- Identified Residues ---")
print("LIM Neighbors:", stored.lim_info)
print("Cys Sites:", stored.cys_info)
python end