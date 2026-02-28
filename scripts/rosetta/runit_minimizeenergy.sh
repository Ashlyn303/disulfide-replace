#!/bin/bash

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

INPUT_DIR="$PROJECT_ROOT/inputs/generated_mutants"
OUTPUT_DIR="$PROJECT_ROOT/results/rosetta_minimizeenergy_results"
SUMMARY_FILE="$PROJECT_ROOT/results/tables/rosetta_minimizeenergy_summary.csv"
ORIGINAL_STRUCT="$PROJECT_ROOT/inputs/5B3N.cif" # Specify either .cif or .pdb path
if [[ "$ORIGINAL_STRUCT" == *.cif ]]; then
    ORIGINAL_FORMAT="CIF"
elif [[ "$ORIGINAL_STRUCT" == *.pdb ]]; then
    ORIGINAL_FORMAT="PDB"
else
    echo "ERROR: Original structure format not recognized (must be .cif or .pdb): $ORIGINAL_STRUCT"
    exit 1
fi
ROSETTA_BIN="$ROSETTA3/bin/minimize.static.linuxgccrelease"
mkdir -p "$OUTPUT_DIR" "$(dirname "$SUMMARY_FILE")"
LOG_DIR="$OUTPUT_DIR"

# --- Group Definitions (Paste EXACTLY from PyMOL Script) ---
PYTHON_GROUPS=$(cat <<EOF
g1_cys = [('A', 'CYS', '24'), ('A', 'CYS', '98')]
g1_lim = [('A', 'LEU', '6'), ('A', 'LEU', '81')]

g2_cys = [('A', 'CYS', '161'), ('A', 'CYS', '230')]
g2_lim = [('A', 'LEU', '142'), ('A', 'MET', '175')]
EOF
)
# ---------------------------------------------------------

# Process Python definitions into Shell variables
eval "$(python3 -c "
import ast
data = {}
exec(\"\"\"$PYTHON_GROUPS\"\"\", {}, data)
def export_group(prefix, lim, cys):
    all_sites = lim + cys
    print(f'{prefix}_COUNT={len(all_sites)}')
    for i, (ch, resn, idx) in enumerate(all_sites, 1):
        print(f'{prefix}_S{i}_CH=\"{ch}\"')
        print(f'{prefix}_S{i}_IDX={idx}')
export_group('G1', data['g1_lim'], data['g1_cys'])
export_group('G2', data['g2_lim'], data['g2_cys'])
")"

echo "filename,group,mutations,total_score" > "$SUMMARY_FILE"

echo "Starting Rosetta batch minimization..."

for pdb in "$INPUT_DIR"/*.pdb; do
    [ -e "$pdb" ] || continue
    basename=$(basename "$pdb" .pdb)
    
    # 1. Identify Group and extract mutations directly from PDB
    if [[ "$basename" == G1* ]]; then
        group="G1"; count=$G1_COUNT
    else
        group="G2"; count=$G2_COUNT
    fi

    muts=""
    for i in $(seq 1 $count); do
        ch_v="${group}_S${i}_CH";   ch="${!ch_v}"
        idx_v="${group}_S${i}_IDX"; idx="${!idx_v}"
        
        # Extract WT from Original Structure
        if [ "$ORIGINAL_FORMAT" == "CIF" ]; then
            # CIF: Field 6 is ResName, 7 is Chain, 9 is ResSeq
            wt_3=$(awk -v ch="$ch" -v idx="$idx" '$1=="ATOM" && $7==ch && $9==idx {print $6; exit}' "$ORIGINAL_STRUCT")
        elif [ "$ORIGINAL_FORMAT" == "PDB" ]; then
            # PDB: Robust fixed-column substr
            wt_3=$(awk -v ch="$ch" -v idx="$idx" '/^ATOM/ && substr($0, 22, 1) == ch && substr($0, 23, 4)+0 == idx {print substr($0, 18, 3); exit}' "$ORIGINAL_STRUCT")
        fi
        wt="${wt_3:0:1}" # Use first letter
        
        # Extract Mutant from PDB (Robust fixed-column substr)
        m=$(awk -v ch="$ch" -v idx="$idx" '/^ATOM/ && substr($0, 22, 1) == ch && substr($0, 23, 4)+0 == idx {print substr($0, 18, 3); exit}' "$pdb")
        
        [ -n "$muts" ] && muts="${muts}-"
        muts="${muts}${wt}${idx}${m}"
    done

    # Clean the PDB for stability
    grep -E "^ATOM|^TER|^END" "$pdb" > "${pdb}.clean"
    
    
    # Run Rosetta and capture output to extract seed
    log_file="$LOG_DIR/${basename}.log"
    $ROSETTA_BIN \
        -s "${pdb}.clean" \
        -score:weights ref2015 \
        -run:min_type lbfgs_armijo \
        -run:min_tolerance 0.0001 \
        -run:constant_seed true \
        -run:jran 11105 \
        -packing:ex1 \
        -packing:ex2 \
        -use_input_sc \
        -flip_HNQ \
        -no_optH false \
        -out:path:all "$OUTPUT_DIR" \
        -out:prefix "${basename}_" \
        -ignore_unrecognized_res \
        -overwrite > "$log_file" 2>&1
    
    # Extract seed from log (specifically looking for seed=VALUE)
    seed=$(grep -i "seed=" "$log_file" | grep -oE "seed=[0-9]+" | head -n 1 | cut -d= -f2)
    
    SCORE_FILE="$OUTPUT_DIR/${basename}_score.sc"
    
    if [ -f "$SCORE_FILE" ]; then
        score=$(grep "SCORE:" "$SCORE_FILE" | tail -n 1 | awk '{print $2}')
    	echo "$basename | Mutations: $muts | Score: $score | Seed: $seed"
        echo "${basename},${group},${muts},${score}" >> "$SUMMARY_FILE"
        # Log file is kept in $OUTPUT_DIR
    else
    	echo "FAILED: $basename | Mutations: $muts | Seed: $seed"
        echo "${basename},${group},${muts},FAILED" >> "$SUMMARY_FILE"
    fi
    
    rm "${pdb}.clean"
done

echo "------------------------------------------------"
echo "Batch processing complete. Summary: $SUMMARY_FILE"