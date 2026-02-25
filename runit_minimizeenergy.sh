#!/bin/bash

# Configuration
INPUT_DIR="generated_mutants"
OUTPUT_DIR="rosetta_minimizeenergy_results"
SUMMARY_FILE="rosetta_minimizeenergy_summary.csv"
ROSETTA_BIN="$ROSETTA3/bin/minimize.static.linuxgccrelease"

# Create output directory and header
mkdir -p "$OUTPUT_DIR"
echo "filename,group,mutations,total_score" > "$SUMMARY_FILE"

echo "Starting Rosetta batch minimization..."

for pdb in "$INPUT_DIR"/*.pdb; do
    [ -e "$pdb" ] || continue
    basename=$(basename "$pdb" .pdb)
    
    # Precise extraction: $6 is the Residue Index, $4 is the Residue Name
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

    # Clean the PDB for stability
    grep -E "^ATOM|^TER|^END" "$pdb" > "${pdb}.clean"
    
    
    # Run Rosetta and capture output to extract seed
    log_file="$OUTPUT_DIR/${basename}.log"
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
    
    # Extract seed from log (highly robust search)
    seed=$(grep -i "seed" "$log_file" | grep -oE "[0-9]{4,}" | head -n 1)
    
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