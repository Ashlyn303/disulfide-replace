reinitialize
python
import os
# Detect project root relative to current working directory
cwd = os.getcwd()
if cwd.endswith("scripts/pymol"):
    project_root = os.path.dirname(os.path.dirname(cwd))
else:
    project_root = cwd

input_pdb = os.path.join(project_root, "inputs/5B3N.cif")
if not os.path.exists(input_pdb):
    print(f"ERROR: Input file {input_pdb} not found. Stopping.")
    quit()

# Load the specific structure
cmd.load(input_pdb, "5B3N")
python end

# Identify Disulfide Pairs and Neighbors
python
from pymol import stored
import itertools

# 1. Get all CYS SG atoms
stored.sg_atoms = []
cmd.iterate("resn CYS and name SG", "stored.sg_atoms.append((chain, resi, index))")

# 2. Find pairs within 2.5A
disulfide_pairs = []
for i in range(len(stored.sg_atoms)):
    for j in range(i + 1, len(stored.sg_atoms)):
        at1 = stored.sg_atoms[i]
        at2 = stored.sg_atoms[j]
        # Use cmd.get_distance for accuracy
        dist = cmd.get_distance(f"index {at1[2]}", f"index {at2[2]}")
        if dist < 2.5:
            disulfide_pairs.append((at1, at2))

# 3. For each pair, find neighbors and build group
groups = []
for idx, (at1, at2) in enumerate(disulfide_pairs, 1):
    cys_residues = []
    # Get full residue info (chain, resn, resi) for the two CYS
    cmd.iterate(f"index {at1[2]} or index {at2[2]}", "cys_residues.append((chain, resn, resi))")
    
    # Selection for this specific pair's neighbors
    pair_sel = f"index {at1[2]} or index {at2[2]}"
    neighbor_sel = f"(byres (resn LEU+ILE+MET and not (name N+CA+C+O)) within 4 of ({pair_sel})) and name CA"
    
    lim_neighbors = []
    cmd.iterate(neighbor_sel, "lim_neighbors.append((chain, resn, resi))")
    
    groups.append({
        'id': f'G{idx}',
        'cys': list(set(cys_residues)), # Unique residues
        'lim': list(set(lim_neighbors))
    })

# 4. Print structured output
print("\n# --- Identified Residues (Copy/Paste this block) ---")
print("Groups = [")
for i, g in enumerate(groups):
    comma = "," if i < len(groups) - 1 else ""
    print("    {")
    print(f"        'id': '{g['id']}',")
    print(f"        'cys': {g['cys']},")
    print(f"        'lim': {g['lim']}")
    print(f"    }}{comma}")
print("]")
print("# --------------------------------------------------")
python end