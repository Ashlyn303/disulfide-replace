# Load the specific structure [cite: 1]
load /Users/szu-hsuan/Library/CloudStorage/Dropbox/python3/anti-HIV_Utag/reimplement_deroo_scFv_disulf_removal_v2/5B3N.cif

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

# Collect (Name, Index) for neighbors and disulfides
cmd.iterate("(near_lim and name CA)", "stored.lim_info.append((resn, resi))")
cmd.iterate("(resn CYS and name CA)", "stored.cys_info.append((resn, resi))")

print("\n--- Identified Residues ---")
print("LIM Neighbors:", stored.lim_info)
print("Cys Sites:", stored.cys_info)
python end