python
import itertools
import os

# Define Groups
g1_cys = [('CYS', '24'), ('CYS', '98')]
g1_lim = [('LEU', '6'), ('LEU', '81')]

g2_cys = [('CYS', '161'), ('CYS', '230')]
g2_lim = [('LEU', '142'), ('MET', '175')]

atom_map = {'CYS': 6, 'ALA': 5, 'VAL': 7, 'LEU': 8, 'ILE': 8, 'MET': 8}
output_dir = "generated_mutants"
if not os.path.exists(output_dir): os.makedirs(output_dir)

def generate_group_mutants(prefix, cys_sites, lim_sites):
    wt_total = sum(atom_map[s[0]] for s in cys_sites + lim_sites)
    all_combos = list(itertools.product(['LEU', 'ILE', 'MET'], ['LEU', 'ILE', 'MET'], ['ALA', 'VAL'], ['ALA', 'VAL']))
    
    count = 0
    for combo in all_combos:
        if sum(atom_map[res] for res in combo) <= wt_total:
            # Re-load clean structure for each mutant
            cmd.load("5B3N.cif", "temp_mutant")
            wizard = cmd.get_wizard()
            
            # Apply 4 mutations (2 LIM, 2 CYS)
            mutations = [
                (lim_sites[0][1], combo[0]), (lim_sites[1][1], combo[1]),
                (cys_sites[0][1], combo[2]), (cys_sites[1][1], combo[3])
            ]
            
            for resi, resn in mutations:
                cmd.wizard("mutagenesis")
                cmd.get_wizard().set_mode(resn)
                cmd.get_wizard().do_select(f"resi {resi}")
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