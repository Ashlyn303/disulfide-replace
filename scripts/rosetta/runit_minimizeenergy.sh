#!/bin/bash

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

INPUT_DIR="$PROJECT_ROOT/inputs/generated_mutants"
OUTPUT_DIR="$PROJECT_ROOT/results/rosetta_minimizeenergy_results"
SUMMARY_FILE="$PROJECT_ROOT/results/tables/rosetta_minimizeenergy_summary.csv"
ORIGINAL_STRUCT="$PROJECT_ROOT/inputs/5B3N.cif" # Specify either .cif or .pdb path
ROSETTA_BIN="$ROSETTA3/bin/minimize.static.linuxgccrelease"
LOG_DIR="$OUTPUT_DIR"

# --- Group Definitions (Paste EXACTLY from PyMOL Script) ---
MUTATION_GROUPS=$(cat <<EOF
Groups = [
    {
        'id': 'G1',
        'cys': [('A', 'CYS', '24'), ('A', 'CYS', '98')],
        'lim': [('A', 'LEU', '6'), ('A', 'LEU', '81')]
    },
    {
        'id': 'G2',
        'cys': [('A', 'CYS', '161'), ('A', 'CYS', '230')],
        'lim': [('A', 'LEU', '142'), ('A', 'MET', '175')]
    }
]
EOF
)
# ---------------------------------------------------------
mkdir -p "$OUTPUT_DIR" "$(dirname "$SUMMARY_FILE")"

if [[ "$ORIGINAL_STRUCT" == *.cif ]]; then
    ORIGINAL_FORMAT="CIF"
elif [[ "$ORIGINAL_STRUCT" == *.pdb ]]; then
    ORIGINAL_FORMAT="PDB"
else
    echo "ERROR: Original structure format not recognized (must be .cif or .pdb): $ORIGINAL_STRUCT"
    exit 1
fi

# Process Python definitions into Shell variables
eval "$(python3 -c "
import ast, sys
data = {}
try:
    exec(\"\"\"$MUTATION_GROUPS\"\"\", {}, data)
    g_key = 'Groups' if 'Groups' in data else 'groups'
    if g_key not in data:
        raise KeyError('Neither \"Groups\" nor \"groups\" found in MUTATION_GROUPS')
    
    groups_list = []
    for G in data[g_key]:
        prefix = G['id']
        groups_list.append(prefix)
        all_sites = G['lim'] + G['cys']
        print(f'{prefix}_COUNT={len(all_sites)}')
        for i, (ch, res_n, idx) in enumerate(all_sites, 1):
            print(f'{prefix}_S{i}_CH=\"{ch}\"')
            print(f'{prefix}_S{i}_IDX={idx}')
    
    # Safe join for bash evaluation
    groups_str = ' '.join(groups_list)
    print(f'ALL_GROUP_IDS=\"{groups_str}\"')
except Exception as e:
    # Use standard formatting to avoid f-string quoting issues
    sys.stderr.write('ERROR in MUTATION_GROUPS parsing: {}\\n'.format(e))
    # Print exit 1 so bash script stops
    print('exit 1')
")"

echo "filename,group,mutations,total_score" > "$SUMMARY_FILE"

echo "Starting Rosetta batch minimization..."

for pdb in "$INPUT_DIR"/*.pdb; do
    [ -e "$pdb" ] || continue
    basename=$(basename "$pdb" .pdb)
    
    # 1. Identify Group and extract mutations directly from PDB
    group=""
    # Use underscore matching to prevent G1 matching G10
    for g_id in $ALL_GROUP_IDS; do
        if [[ "$basename" == ${g_id}_* ]]; then
            group="$g_id"; break
        fi
    done

    if [ -z "$group" ]; then
        echo "WARNING: Skipping $basename - does not match any group prefix with underscore (e.g., G1_). Current ALL_GROUP_IDS: '$ALL_GROUP_IDS'" >&2
        continue
    fi

    count_var="${group}_COUNT"
    count=${!count_var}

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
        m_3=$(awk -v ch="$ch" -v idx="$idx" '/^ATOM/ && substr($0, 22, 1) == ch && substr($0, 23, 4)+0 == idx {print substr($0, 18, 3); exit}' "$pdb")
        m="${m_3:0:1}" # Use first letter
        
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