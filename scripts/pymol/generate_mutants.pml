reinitialize
python
import itertools
import os

# Detect project root relative to current working directory
cwd = os.getcwd()
if cwd.endswith("scripts/pymol"):
    project_root = os.path.dirname(os.path.dirname(cwd))
else:
    project_root = cwd

# Define Groups
g1_cys = [('A', 'CYS', '24'), ('A', 'CYS', '98')]
g1_lim = [('A', 'LEU', '6'), ('A', 'LEU', '81')]

g2_cys = [('A', 'CYS', '161'), ('A', 'CYS', '230')]
g2_lim = [('A', 'LEU', '142'), ('A', 'MET', '175')]

atom_map = {'CYS': 6, 'ALA': 5, 'VAL': 7, 'LEU': 8, 'ILE': 8, 'MET': 8}
output_dir = os.path.join(project_root, "inputs/generated_mutants")
if not os.path.exists(output_dir): os.makedirs(output_dir)

def generate_group_mutants(prefix, cys_sites, lim_sites):
    wt_total = sum(atom_map[s[1]] for s in cys_sites + lim_sites)
    
    # LIM sites mutate to L/I/M, CYS sites mutate to A/V
    lim_options = [['LEU', 'ILE', 'MET']] * len(lim_sites)
    cys_options = [['ALA', 'VAL']] * len(cys_sites)
    all_combos = list(itertools.product(*(lim_options + cys_options)))
    
    count = 0
    for combo in all_combos:
        if sum(atom_map[res] for res in combo) <= wt_total:
            # Re-load clean structure for each mutant using detected root
            input_pdb = os.path.join(project_root, "inputs/5B3N.cif")
            if not os.path.exists(input_pdb):
                print(f"ERROR: Input file {input_pdb} not found. Stopping.")
                return
            cmd.delete("temp_mutant")
            cmd.load(input_pdb, "temp_mutant")
            wizard = cmd.get_wizard()
            
            # Apply mutations to all identified sites
            all_sites = lim_sites + cys_sites
            for (chain, resn, resi), new_resn in zip(all_sites, combo):
                cmd.wizard("mutagenesis")
                cmd.get_wizard().set_mode(new_resn)
                cmd.get_wizard().do_select(f"/temp_mutant//{chain}/{resi}")
                cmd.get_wizard().apply()
            
            filename = f"{output_dir}/{prefix}_mutant_{count:03d}.pdb"
            cmd.save(filename, "temp_mutant")
            cmd.delete("temp_mutant")
            count += 1
    print(f"Generated {count} mutants for {prefix}")

# Run for both groups independently
generate_group_mutants("G1", g1_cys, g1_lim)
generate_group_mutants("G2", g2_cys, g2_lim)
python end