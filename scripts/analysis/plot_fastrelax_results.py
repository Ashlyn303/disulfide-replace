#%%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import os
import re

# Configuration: Filenames relative to the script's directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# New: Handle nested directory structure in scripts/analysis/
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
INPUT_FILE = os.path.join(PROJECT_ROOT, "results", "tables", "rosetta_fastrelax_summary.csv")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "results", "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)

AA_MAP = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def shorten_mutations(mut_str):
    """Convert L6LEU-L81LEU to L6L-L81L formatting."""
    parts = str(mut_str).split('-')
    short_parts = []
    for p in parts:
        # p is like 'L6LEU'
        # We find where the 3-letter code starts
        for aa3 in AA_MAP:
            if p.endswith(aa3):
                prefix = p[:-3]
                short_parts.append(f"{prefix}{AA_MAP[aa3]}")
                break
        else:
            short_parts.append(p) # fallback
    return "-".join(short_parts)

def extract_mutant_id(filename):
    """Clean 'G1_mutant_000_G1_mutant_000.pdb_0001' -> 'G1_mutant_000'."""
    # Matches G1_mutant_000 or Group1_mutant_000 etc.
    match = re.search(r'(Group\d+|G\d+)_mutant_\d+', str(filename))
    if match:
        return match.group(0)
    return str(filename)

def set_plot_style(small_size=18, medium_size=22, large_size=26, 
                   figsize=(12, 8), linewidth=3, tick_width=2, tick_length=8,
                   labelpad=10):
    """
    Set comprehensive matplotlib styling. Call this function at the start 
    of any script to apply consistent styling.
    """
    plt.rc('font', size=small_size)
    plt.rc('axes', titlesize=medium_size)
    plt.rc('axes', labelsize=medium_size)
    plt.rc('xtick', labelsize=medium_size)
    plt.rc('ytick', labelsize=medium_size)
    plt.rc('legend', fontsize=small_size)
    plt.rc('figure', titlesize=small_size)
    plt.rc('lines', linewidth=linewidth)
    plt.rc('lines', markersize=10)
    plt.rc('axes', linewidth=2)
    plt.rc('axes', edgecolor='black')  # Black borders
    
    # Figure settings
    mpl.rcParams['figure.figsize'] = figsize
    mpl.rcParams['figure.dpi'] = 100
    mpl.rcParams['savefig.dpi'] = 300
    
    # Grid settings
    mpl.rcParams['axes.grid'] = True
    mpl.rcParams['grid.alpha'] = 0.3
    
    # Tick settings - thickness, length, and direction
    mpl.rcParams['xtick.major.width'] = tick_width
    mpl.rcParams['ytick.major.width'] = tick_width
    mpl.rcParams['xtick.minor.width'] = tick_width * 0.6
    mpl.rcParams['ytick.minor.width'] = tick_width * 0.6
    mpl.rcParams['xtick.major.size'] = tick_length
    mpl.rcParams['ytick.major.size'] = tick_length
    mpl.rcParams['xtick.minor.size'] = tick_length * 0.6
    mpl.rcParams['ytick.minor.size'] = tick_length * 0.6
    mpl.rcParams['xtick.color'] = 'black'  # Black tick marks
    mpl.rcParams['ytick.color'] = 'black'  # Black tick marks
    
    # Tick direction: 'in', 'out', or 'inout'
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    
    # Show ticks on all sides for complete box
    mpl.rcParams['xtick.top'] = False
    mpl.rcParams['xtick.bottom'] = True
    mpl.rcParams['ytick.left'] = True
    mpl.rcParams['ytick.right'] = False
    
    # Axis label padding - creates space between labels and axes
    mpl.rcParams['axes.labelpad'] = labelpad

def main():
    set_plot_style() # Apply custom styling
    
    if not os.path.exists(INPUT_FILE):
        print(f"Error: {INPUT_FILE} not found.")
        return

    # Load data
    df = pd.read_csv(INPUT_FILE)
    
    # Filter out FAILED rows
    df = df[df['avg_score'] != 'FAILED'].copy()
    
    # Dynamically identify all replicate columns (rep1, rep2, ..., repN)
    rep_cols = [col for col in df.columns if col.startswith('rep')]
    for col in rep_cols + ['avg_score', 'std_score']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Get unique groups to generate separate figures
    unique_groups = sorted(df['group'].unique())
    print(f"Found groups: {unique_groups}")

    for group in unique_groups:
        print(f"Generating plot for group: {group}...")
        
        # Filter data for this group
        group_df = df[df['group'] == group].copy()
        
        # Create refined labels: "L6L-L81L... \n (G1_mutant_000)"
        group_df['plot_label'] = group_df.apply(
            lambda x: f"{shorten_mutations(x['mutations'])}\n({extract_mutant_id(x['filename'])})", axis=1
        )
        
        # Sort by avg_score (most negative first)
        group_df = group_df.sort_values('avg_score', ascending=True)

        # Melt the dataframe for plotting individual replicates
        group_df_melted = group_df.melt(
            id_vars=['plot_label', 'group', 'avg_score'], 
            value_vars=rep_cols, 
            var_name='replicate', 
            value_name='score'
        )

        # Plotting
        plt.figure(figsize=(10, 8)) # Narrower figure for a more compact x-axis
        
        # Enforce the sorting order explicitly
        sort_order = group_df['plot_label'].tolist()

        # 1. Plot the Box Plot (Showing distribution: median, quartiles, range)
        sns.boxplot(
            data=group_df_melted,
            x='plot_label',
            y='score',
            order=sort_order,
            palette='viridis',
            width=0.4, # More compact boxes
            showfliers=False,
            showmeans=True,
            meanline=True,
            meanprops={'color': 'red', 'ls': '--', 'lw': 1.5}, # Red dashed line for mean
            medianprops={'color': 'blue', 'lw': 2},          # Blue solid line for median
            boxprops=dict(alpha=0.6)
        )

        # 2. Overlay individual data points (Strip plot)
        sns.stripplot(
            data=group_df_melted,
            x='plot_label',
            y='score',
            order=sort_order,
            color='black',
            size=3, # Slightly smaller points for thinner boxes
            jitter=True,
            alpha=0.4
        )

        # 3. Add terminal output for numerical mean and median
        print(f"\n--- Statistics for Group {group} ---")
        # We still use the original 'mutations' for the print table for clarity if needed, 
        # or we can use the new plot_label. Let's stick to the combined label for consistency.
        for i, row in group_df.iterrows():
            mut_label = row['plot_label'].replace('\n', ' ')
            avg = row['avg_score']
            # Re-calculate median from melted if needed, or just use avg for simple table
            mut_data = group_df_melted[group_df_melted['plot_label'] == row['plot_label']]['score'].dropna()
            med = mut_data.median()
            print(f"{mut_label:50} | Mean: {avg:8.2f} | Median: {med:8.2f}")

        # Determine output filename for this group
        group_output = os.path.join(OUTPUT_DIR, f"fastrelax_results_{group}_plot.png")

        # plt.title(f'Rosetta FastRelax Results - Group {group} (Sorted by Mean Score)')
        plt.xlabel('Mutantations')
        plt.ylabel('REU')
        plt.ylim(-780, -755)
        
        # Add a custom legend for mean/median lines
        from matplotlib.lines import Line2D
        custom_lines = [Line2D([0], [0], color='red', lw=2, ls='--'),
                        Line2D([0], [0], color='blue', lw=2)]
        plt.legend(custom_lines, ['Mean', 'Median'], loc='lower right')

        # Rotate x-axis labels for better readability
        plt.xticks(rotation=45, ha='right', fontsize=14)
        
        plt.tight_layout()
        plt.savefig(group_output, dpi=300)
        print(f"Plot saved to {group_output}")
        plt.show() # Display the plot interactively
        plt.close() # Close figure to free memory

if __name__ == "__main__":
    main()

# %%
