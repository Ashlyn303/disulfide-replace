#!/bin/bash

# Configuration
INPUT_DIR="generated_mutants"
OUTPUT_DIR="rosetta_fastrelax_results"
SUMMARY_FILE="rosetta_fastrelax_summary.csv"
ROSETTA_BIN="$ROSETTA3/bin/relax.static.linuxgccrelease"
REPLICATES=3

# Create output directory and header
mkdir -p "$OUTPUT_DIR"
echo "filename,group,mutations,rep1,rep2,rep3,avg_score" > "$SUMMARY_FILE"

echo "Starting Rosetta FastRelax (3 replicates per mutant)..."

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
    
    scores=()
    for i in $(seq 1 $REPLICATES); do
        echo "Processing: $basename (Rep $i) | Mutations: $muts"
        
        # Run FastRelax
        # -relax:fast: standard 5-cycle relax
        # -nstruct 1: we handle the loop ourselves for specific naming
        $ROSETTA_BIN \
            -s "${pdb}.clean" \
            -relax:fast \
            -out:path:all "$OUTPUT_DIR" \
            -out:prefix "${basename}_rep${i}_" \
            -overwrite > /dev/null 2>&1
        
        SCORE_FILE="$OUTPUT_DIR/${basename}_rep${i}_score.sc"
        if [ -f "$SCORE_FILE" ]; then
            val=$(grep "SCORE:" "$SCORE_FILE" | tail -n 1 | awk '{print $2}')
            scores+=($val)
        else
            scores+=("FAILED")
        fi
    done

    # 3. Calculate Average (sum divided by number of successful replicates)
    avg=$(python3 -c "import sys; s=[float(x) for x in sys.argv[1:] if x!='FAILED']; print(round(sum(s)/len(s), 3) if s else 'FAILED')" "${scores[@]}")
    
    # 4. Append to CSV
    echo "${basename},${group},${muts},${scores[0]},${scores[1]},${scores[2]},$avg" >> "$SUMMARY_FILE"
    
    rm "${pdb}.clean"
done

echo "------------------------------------------------"
echo "Relaxation complete. Summary: $SUMMARY_FILE"
