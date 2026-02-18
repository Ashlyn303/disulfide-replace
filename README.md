# Disulfide Replace

This project provides a workflow for identifying disulfide bonds in protein structures, finding their neighboring residues, and performing combinatorial mutations to replace the disulfide bonds with hydrophobic packing (Ala/Val) and adjusting neighbors (Leu/Ile/Met). The resulting mutants are then evaluated using Rosetta.

## Technical Rationale

The goal is to search within maximal heavy atom conserving and hydrophobic combinations to replace existing disulfide bonds.

1.  **Neighbor Identification**: Identify disulfide neighbor Leu, Ile, Met residues that have a sidechain atom (not N, CA, C, O) within 4Å of the disulfide sidechains (CB, SG, CB, SG).
2.  **Combinatorial Search**: Construct a search where we consider **Ala and Val** at each former Cys site, and consider **Leu, Ile, and Met** at each neighboring site that is already Leu, Ile, or Met.
    -   Example: A local cluster of `[Leu, Cys, Met, Cys]` becomes `[LIM][AV][LIM][AV]`.
3.  **Heavy Atom Conservation**: Within this set, only combinations that **maintain or reduce the number of heavy atoms** are legal (e.g., no double Valine mutants).

## Workflow Overview

The workflow consists of four main steps:

1.  **Structure Preparation**: The input `5B3N.pdb` is downloaded directly from [RCSB PDB](https://www.rcsb.org/structure/5B3N).
2.  **Neighbor Identification**: Finding residue positions near the disulfide sites.
3.  **Combinatorial Mutagenesis**: Generating all valid mutant combinations using PyMOL.
4.  **Rosetta Evaluation**: Scoring and relaxing the generated mutants.

## Prerequisites

-   **PyMOL**: Required for running `.pml` scripts.
-   **Rosetta**: Required for energy minimization and FastRelax calculations. Ensure `$ROSETTA3` environment variable is set.
-   **Python 3**: For data processing within PyMOL and shell scripts.

## Step-by-Step Instructions

### 1. Identify Neighbors
Run `search_within_4A.pml` in PyMOL to identify Leu, Ile, and Met residues within 4Å of the disulfide sidechains (CB, SG).

```bash
pymol -qc search_within_4A.pml
```

This identifies two groups in `5B3N`:
- **Group 1**: Cys 24, 98 | Neighbors Leu 6, 81
- **Group 2**: Cys 161, 230 | Neighbors Leu 142, Met 175

### 2. Generate Mutants
Run `generate_mutants.pml` to perform combinatorial mutagenesis based on the rationale above.

```bash
pymol -qc generate_mutants.pml
```

This will create a `generated_mutants` folder containing the mutated `.pdb` files.

### 3. Rosetta Evaluation
Run the provided shell scripts to evaluate the mutants using Rosetta.

#### FastRelax
Performs a standard FastRelax with 3 replicates per mutant and generates a summary CSV.

```bash
./runit_fastrelax.sh
```

#### Energy Minimization
Performs simple energy minimization on the mutants.

```bash
./runit_minimizeenergy.sh
```

## Outputs

-   **`generated_mutants/`**: Folder containing mutated PDB structures.
-   **`rosetta_fastrelax_results/`**: Rosetta FastRelax output files.
-   **`rosetta_minimizeenergy_results/`**: Rosetta minimization output files.
-   **`rosetta_fastrelax_summary.csv`**: Aggregated scores and mutation details for FastRelax.
-   **`rosetta_minimizeenergy_summary.csv`**: Aggregated scores for minimization.

## Software Citations

### Rosetta
"Computational modeling and design were performed using Rosetta version 2023.45 (release a6d9ba8). Leman, J. K., et al. (2020). 'Macromolecular modeling and design in Rosetta: recent methods and frameworks.' Nature Methods, 17(7), 665-680."

### PyMOL
"Molecular visualizations were generated using The PyMOL Molecular Graphics System, Version 3.1.3, Schrödinger, LLC."
