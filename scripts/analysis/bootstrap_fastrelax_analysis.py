#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

# Configuration
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "tables")
FIGURES_DIR = os.path.join(PROJECT_ROOT, "results", "figures")
INPUT_FILES = [
    os.path.join(TABLES_DIR, "rosetta_fastrelax_summary-1.csv"),
    os.path.join(TABLES_DIR, "rosetta_fastrelax_summary-2.csv"),
    os.path.join(TABLES_DIR, "rosetta_fastrelax_summary-3.csv")
]
OUTPUT_CSV = os.path.join(TABLES_DIR, "rosetta_fastrelax_bootstrap_results.csv")
OUTPUT_PLOT = os.path.join(FIGURES_DIR, "fastrelax_bootstrap_comparison.png")

os.makedirs(FIGURES_DIR, exist_ok=True)

def bootstrap_stats(data, n_bootstrap=1000, ci=0.95):
    """Perform bootstrap resampling to estimate mean, std, and CI."""
    data = np.array(data)
    data = data[~np.isnan(data)] # Remove NaNs
    if len(data) == 0:
        return np.nan, np.nan, np.nan, np.nan
    
    boot_means = []
    for _ in range(n_bootstrap):
        resample = np.random.choice(data, size=len(data), replace=True)
        boot_means.append(np.mean(resample))
    
    boot_mean = np.mean(boot_means)
    boot_std = np.std(boot_means)
    
    alpha = 1.0 - ci
    lower = np.percentile(boot_means, (alpha/2.0) * 100)
    upper = np.percentile(boot_means, (1.0 - alpha/2.0) * 100)
    return boot_mean, boot_std, lower, upper

def main():
    all_data = []
    for f in INPUT_FILES:
        if os.path.exists(f):
            print(f"Loading {f}...")
            df = pd.read_csv(f)
            all_data.append(df)
        else:
            print(f"Warning: {f} not found.")
    
    if not all_data:
        print("Error: No input files found.")
        return

    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Identify unique mutants by filename and mutations
    # Filenames might have slight variations if they were from different runs but should be the same mutants
    # Let's clean the filename to be safe
    def clean_fn(fn):
         match = re.search(r'(Group\d+|G\d+)_mutant_\d+', str(fn))
         return match.group(0) if match else str(fn)
    
    combined_df['clean_name'] = combined_df['filename'].apply(clean_fn)
    
    # Identify replicate columns
    rep_cols = [c for c in combined_df.columns if c.startswith('rep')]
    
    results = []
    
    # Group by clean name and mutations
    grouped = combined_df.groupby(['clean_name', 'mutations', 'group'])
    
    for (name, muts, group), group_data in grouped:
        # Collect all replicates for this mutant
        all_reps = group_data[rep_cols].values.flatten()
        all_reps = pd.to_numeric(all_reps, errors='coerce')
        
        b_mean, b_std, b_lower, b_upper = bootstrap_stats(all_reps)
        orig_mean = np.nanmean(all_reps)
        orig_std = np.nanstd(all_reps)
        
        results.append({
            'name': name,
            'mutations': muts,
            'group': group,
            'original_mean': orig_mean,
            'original_std': orig_std,
            'bootstrap_mean': b_mean,
            'bootstrap_std': b_std,
            'ci_lower': b_lower,
            'ci_upper': b_upper,
            'n_samples': len(all_reps[~np.isnan(all_reps)])
        })
    
    res_df = pd.DataFrame(results)
    res_df.to_csv(OUTPUT_CSV, index=False)
    print(f"Bootstrap results saved to {OUTPUT_CSV}")
    
    # Print results for each group
    for group_name in sorted(res_df['group'].unique()):
        print(f"\n" + "="*80)
        print(f" BOOTSTRAP RESULTS FOR GROUP: {group_name}")
        print("="*80)
        print(f"{'Mutant':<25} | {'Mean ± STD':>15} | {'95% CI Lower':>15} | {'95% CI Upper':>15} | {'N':>5}")
        print("-" * 80)
        
        group_results = res_df[res_df['group'] == group_name].sort_values('bootstrap_mean')
        for _, row in group_results.iterrows():
            mean_std = f"{row['original_mean']:8.2f} ± {row['original_std']:5.2f}"
            print(f"{row['name']:<25} | {mean_std:>15} | {row['ci_lower']:15.2f} | {row['ci_upper']:15.2f} | {int(row['n_samples']):>5}")

    # Plotting (Disabled per user request)
    # plt.figure(figsize=(12, 8))
    # sns.set_style("whitegrid")
    # ... (code truncated for brevity)
    # plt.savefig(OUTPUT_PLOT, dpi=300)
    # print(f"Comparison plot saved to {OUTPUT_PLOT}")

if __name__ == "__main__":
    main()
