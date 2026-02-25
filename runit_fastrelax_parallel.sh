#!/bin/bash

# Configuration
INPUT_DIR="rosetta_minimizeenergy_results_top5"
#INPUT_DIR="tmp"
OUTPUT_DIR="rosetta_fastrelax_results"
SUMMARY_FILE="rosetta_fastrelax_summary.csv"
ROSETTA_BIN="$ROSETTA3/bin/relax.static.linuxgccrelease"
REPLICATES=20
PARALLEL_JOBS=20

# Create output directory and header
mkdir -p "$OUTPUT_DIR"
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
    if [[ "$basename" == G1* ]]; then
        group="G1"
        m1=$(awk '$6 == 6 {print $4; exit}' "$pdb")
        m2=$(awk '$6 == 81 {print $4; exit}' "$pdb")
        m3=$(awk '$6 == 24 {print $4; exit}' "$pdb")
        m4=$(awk '$6 == 98 {print $4; exit}' "$pdb")
        muts="L6${m1}-L81${m2}-C24${m3}-C98${m4}"
    else
        group="G2"
        m1=$(awk '$6 == 142 {print $4; exit}' "$pdb")
        m2=$(awk '$6 == 175 {print $4; exit}' "$pdb")
        m3=$(awk '$6 == 161 {print $4; exit}' "$pdb")
        m4=$(awk '$6 == 230 {print $4; exit}' "$pdb")
        muts="L142${m1}-M175${m2}-C161${m3}-C230${m4}"
    fi

    # 2. Clean PDB
    grep -E "^ATOM|^TER|^END" "$pdb" > "${pdb}.clean"

    tmpdir=$(mktemp -d)

    # 3. Run FastRelax replicates in parallel
    for i in $(seq 1 $REPLICATES); do
        (
            log_file="$OUTPUT_DIR/${basename}_rep${i}.log"
            $ROSETTA_BIN \
                -s "${pdb}.clean" \
                -relax:fast \
                -run:constant_seed false \
                -out:path:all "$OUTPUT_DIR" \
                -out:prefix "${basename}_rep${i}_" \
                -overwrite > "$log_file" 2>&1

            # Extract seed (highly robust search)
            seed=$(grep -i "seed" "$log_file" | grep -oE "[0-9]{4,}" | head -n 1)
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
