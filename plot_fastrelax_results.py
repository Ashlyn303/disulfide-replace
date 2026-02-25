#%%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Configuration: Filenames relative to the script's directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_FILE = os.path.join(SCRIPT_DIR, "rosetta_fastrelax_summary.csv")
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "fastrelax_results_plot.png")

def main():
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
        
        # Sort by avg_score (most negative first)
        group_df = group_df.sort_values('avg_score', ascending=True)

        # Melt the dataframe for plotting individual replicates
        group_df_melted = group_df.melt(
            id_vars=['mutations', 'group', 'avg_score'], 
            value_vars=rep_cols, 
            var_name='replicate', 
            value_name='score'
        )

        # Plotting
        plt.figure(figsize=(14, 8))
        sns.set_theme(style="whitegrid")
        
        # Enforce the sorting order explicitly
        sort_order = group_df['mutations'].tolist()

        # 1. Plot the Box Plot (Showing distribution: median, quartiles, range)
        sns.boxplot(
            data=group_df_melted,
            x='mutations',
            y='score',
            order=sort_order,
            palette='viridis',
            width=0.6,
            showfliers=False,
            boxprops=dict(alpha=0.6)
        )

        # Determine output filename for this group
        group_output = os.path.join(SCRIPT_DIR, f"fastrelax_results_{group}_plot.png")

        plt.title(f'Rosetta FastRelax Results - Group {group} (Sorted by Mean Score)', fontsize=16)
        plt.xlabel('Mutant (Mutations)', fontsize=12)
        plt.ylabel('Rosetta Total Score (ref2015)', fontsize=12)
        
        # Rotate x-axis labels for better readability
        plt.xticks(rotation=45, ha='right')
        
        plt.tight_layout()
        plt.savefig(group_output, dpi=300)
        print(f"Plot saved to {group_output}")
        plt.show() # Display the plot interactively
        plt.close() # Close figure to free memory

if __name__ == "__main__":
    main()

# %%
