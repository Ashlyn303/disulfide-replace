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
    
    # Ensure scores are numeric
    rep_cols = [f'rep{i}' for i in range(1, 21)]
    for col in rep_cols + ['avg_score', 'std_score']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Sort by avg_score (most negative first)
    df = df.sort_values('avg_score', ascending=True)

    # Melt the dataframe for plotting individual replicates
    # This transforms the data into a 'long' format suitable for sns.stripplot
    df_melted = df.melt(
        id_vars=['mutations', 'group', 'avg_score'], 
        value_vars=rep_cols, 
        var_name='replicate', 
        value_name='score'
    )

    # Plotting
    plt.figure(figsize=(14, 8))
    sns.set_theme(style="whitegrid")
    
    # 1. Plot the Box Plot (Showing distribution: median, quartiles, range)
    # This naturally handles the 'distribution' without drawing bars from zero
    sns.boxplot(
        data=df_melted,
        x='mutations',
        y='score',
        hue='group',
        palette='viridis',
        dodge=False,
        width=0.6,
        showfliers=False,  # Fliers are shown by the stripplot
        boxprops=dict(alpha=0.6)
    )

    # 2. Overlay individual points (Stripplot)
    # This shows every single one of the 20 replicates
    # sns.stripplot(
    #     data=df_melted,
    #     x='mutations',
    #     y='score',
    #     hue='group',
    #     dodge=False,
    #     jitter=0.2,
    #     size=4,
    #     alpha=0.5,
    #     palette='viridis',
    #     legend=False  # Avoid duplicate legend
    # )

    # plt.title('Rosetta FastRelax Results (Score Distributions)', fontsize=16)
    # plt.xlabel('Mutant (Mutations)', fontsize=12)
    # plt.ylabel('Rosetta Total Score (ref2015)', fontsize=12)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    
    # Legend handling
    plt.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_FILE, dpi=300)
    print(f"Plot saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
