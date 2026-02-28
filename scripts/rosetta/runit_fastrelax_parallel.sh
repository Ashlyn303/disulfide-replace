#!/bin/bash

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

INPUT_DIR="$PROJECT_ROOT/results/rosetta_minimizeenergy_results_top5"
OUTPUT_DIR="$PROJECT_ROOT/results/rosetta_fastrelax_results"
SUMMARY_FILE="$PROJECT_ROOT/results/tables/rosetta_fastrelax_summary.csv"
ORIGINAL_STRUCT="$PROJECT_ROOT/inputs/5B3N.cif" # Specify either .cif or .pdb path
ROSETTA_BIN="$ROSETTA3/bin/relax.static.linuxgccrelease"
REPLICATES=20
PARALLEL_JOBS=20

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

# Create output directory and header
mkdir -p "$OUTPUT_DIR" "$(dirname "$SUMMARY_FILE")"
LOG_DIR="$OUTPUT_DIR"
header="filename,group,mutations"
for i in $(seq 1 $REPLICATES); do
    header="${header},rep${i}"
done
header="${header},avg_score,std_score"
echo "$header" > "$SUMMARY_FILE"

echo "Starting Rosetta FastRelax (${REPLICATES} replicates, ${PARALLEL_JOBS} jobs in parallel for each mutant in $INPUT_DIR)..."

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

    # 2. Clean PDB
    grep -E "^ATOM|^TER|^END" "$pdb" > "${pdb}.clean"

    tmpdir=$(mktemp -d)

    # 3. Run FastRelax replicates in parallel
    for i in $(seq 1 $REPLICATES); do
        (
            log_file="$LOG_DIR/${basename}_rep${i}.log"
            $ROSETTA_BIN \
                -s "${pdb}.clean" \
                -relax:fast \
                -score:weights ref2015 \
                -run:constant_seed false \
                -out:path:all "$OUTPUT_DIR" \
                -out:prefix "${basename}_rep${i}_" \
                -overwrite > "$log_file" 2>&1

            # Extract seed from log (specifically looking for seed=VALUE)
            sync
            seed=$(grep -ai "seed=" "$log_file" | grep -oE "seed=[-]?[0-9]+" | head -n 1 | cut -d= -f2)
            echo "Processing: $basename (Rep $i) | Mutations: $muts | Seed: $seed"

            SCORE_FILE="$OUTPUT_DIR/${basename}_rep${i}_score.sc"
            if [ -f "$SCORE_FILE" ]; then
                val=$(grep "SCORE:" "$SCORE_FILE" | tail -n 1 | awk '{print $2}')
                echo "$val" > "$tmpdir/rep${i}.score"
            else
                echo "FAILED" > "$tmpdir/rep${i}.score"
            fi
        ) &

        # Throttle to PARALLEL_JOBS
        if (( $(jobs -pr | wc -l) >= PARALLEL_JOBS )); then
            wait -n
        fi
    done

    wait

    # 4. Collect scores in order
    scores=()
    for i in $(seq 1 $REPLICATES); do
        if [ -f "$tmpdir/rep${i}.score" ]; then
            scores+=("$(cat "$tmpdir/rep${i}.score")")
        else
            scores+=("FAILED")
        fi
    done

    # 5. Calculate Average and Std Dev using Python
    stats=$(python3 -c "
import sys, statistics
data = [float(x) for x in sys.argv[1:] if x != 'FAILED']
if not data:
    print('FAILED,FAILED')
else:
    avg = round(statistics.mean(data), 3)
    std = round(statistics.stdev(data), 3) if len(data) > 1 else 0.0
    print(f'{avg},{std}')
" "${scores[@]}")

    # 6. Append to CSV
    score_line=$(IFS=,; echo "${scores[*]}")
    echo "${basename},${group},${muts},${score_line},${stats}" >> "$SUMMARY_FILE"

    rm -rf "$tmpdir" "${pdb}.clean"
done

echo "------------------------------------------------"
echo "Relaxation complete. Summary: $SUMMARY_FILE"
