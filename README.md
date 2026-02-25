# Disulfide Replace: Generative Design of Peptides using Rosetta

**Elevator Pitch**: This repository contains the computational pipeline for the **mBG17 project**, focusing on the generative design of peptides using Rosetta. It automates the transition from sequence generation to high-resolution energy minimization and statistical stability ranking, enabling the replacement of disulfide bonds with optimized hydrophobic packing.

## Quick Start
Since Rosetta environments are complex, we provide streamlined entry points. Assuming Rosetta is installed and `$ROSETTA3` is set:

```bash
# 1. Identify disulfide neighbors (Group 1 & 2)
pymol -qc scripts/pymol/search_within_4A.pml

# 2. Generate combinatorial hydrophobic mutants
pymol -qc scripts/pymol/generate_mutants.pml

# 3. High-resolution scoring and relaxation (parallel)
./scripts/rosetta/runit_fastrelax_parallel.sh

# 4. Visualize results
python scripts/analysis/plot_fastrelax_results.py
```

## Reproducibility Note
To ensure **100% reproducibility** of the results presented in our research, all energy minimization simulations were performed with a constant seed:
> **Magic Seed: `-run:jran 11105`**

This prevents stochastic variance from affecting the final mutant rankings during the initial optimization phase. Note that FastRelax replicates use random seeds to ensure diverse sampling and robust convergence.

## Repository Organization

The project is structured to separate inputs, logic, and results for a professional scientific workflow:

- **`inputs/`**: Source protein structures (e.g., `5B3N.cif`) and **generated mutants**.
- **`scripts/`** (Pipeline Logic):
    - `pymol/`: Identification and mutagenesis scripts (`.pml`).
    - `rosetta/`: Production Rosetta execution scripts (`.sh`).
    - `analysis/`: Plotting and statistics scripts (`.py`).
- **`results/`**: 
    - `figures/`: Professional plots and visualizations.
    - `tables/`: Summary CSV files with aggregated rankings.
    - `rosetta_fastrelax_results/`: Relaxation output and **diagnostic logs**.
    - `rosetta_minimizeenergy_results/`: Minimization output and **diagnostic logs**.
- **`metadata`**:
    - `DATA.md`: Detailed data dictionary for CSV columns and REU units.
    - `requirements.txt`: Python environment specifications.

## Scientific Workflow

1.  **Neighbor Identification**: Identify disulfide neighbor Leu, Ile, Met residues within 4Å of the disulfide sidechains.
2.  **Combinatorial Search**: Construct a search considering **Ala/Val** at Cys sites and **Leu/Ile/Met** at neighboring sites.
3.  **Hydrophobic Optimization**: Combinations are filtered to maintain or reduce total heavy atom count relative to wildtype.
4.  **Ensemble Evaluation**: Mutants are ranked by the mean score across an ensemble of 20 FastRelax replicates.

## Software Citations

### Rosetta
"Computational modeling and design were performed using Rosetta version 2023.45 (release a6d9ba8). Leman, J. K., et al. (2020)."

### PyMOL
"Molecular visualizations were generated using The PyMOL Molecular Graphics System, Version 3.1.3, Schrödinger, LLC."
