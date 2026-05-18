#!/usr/bin/env bash
# =============================================================================
# PICRUSt2: Installation and full pipeline run
# =============================================================================
# Input  : QIIME2 artifacts (table.qza, rep-seqs.qza)
# Output : picrust2_output/ (MetaCyc pathways, EC numbers, KO terms)
#          stratified = per-ASV contributions included
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CONDA_BASE="$HOME/miniconda3"
PICRUST2_ENV="$HOME/Desktop/SW/picrust2-env"
THREADS=4          # adjust to available CPUs
INPUT_DIR="picrust2_input"
OUTPUT_DIR="picrust2_output"

# ---------------------------------------------------------------------------
# STEP 1 — Install PICRUSt2 into ~/Desktop/SW/picrust2-env
# Skip if already installed
# ---------------------------------------------------------------------------
source "$CONDA_BASE/etc/profile.d/conda.sh"

if [ -d "$PICRUST2_ENV" ]; then
    echo "[INFO] PICRUSt2 environment already exists at $PICRUST2_ENV — skipping install."
else
    echo "[INFO] Creating PICRUSt2 conda environment at $PICRUST2_ENV ..."
    echo "[INFO] This downloads ~2 GB and may take 20-40 minutes with conda."

    # Install mamba into base first for much faster dependency solving
    conda install -n base -c conda-forge mamba -y --quiet 2>/dev/null || true

    if command -v mamba &>/dev/null; then
        SOLVER=mamba
    else
        SOLVER=conda
    fi

    $SOLVER create \
        --prefix "$PICRUST2_ENV" \
        -c conda-forge -c bioconda \
        picrust2=2.6.3 \
        -y

    echo "[INFO] PICRUSt2 installed successfully."
fi

conda activate "$PICRUST2_ENV"
echo "[INFO] PICRUSt2 version: $(picrust2_pipeline.py --version 2>&1 | head -1)"

# ---------------------------------------------------------------------------
# STEP 2 — Extract QIIME2 artifacts (qza = zip files)
# ---------------------------------------------------------------------------
mkdir -p "$INPUT_DIR"

echo "[INFO] Extracting feature table from table.qza ..."
unzip -p table.qza "*/data/feature-table.biom" > "$INPUT_DIR/feature-table.biom"

echo "[INFO] Extracting ASV sequences from rep-seqs.qza ..."
unzip -p rep-seqs.qza "*/data/dna-sequences.fasta" > "$INPUT_DIR/dna-sequences.fasta"

N_ASV=$(grep -c "^>" "$INPUT_DIR/dna-sequences.fasta")
echo "[INFO] ${N_ASV} ASVs extracted."

# ---------------------------------------------------------------------------
# STEP 3 — Run PICRUSt2 full pipeline
# --stratified  : include per-ASV contributions (needed for family-level analysis)
# --in_traits   : predict EC, KO, and pathways (default)
# ---------------------------------------------------------------------------
if [ -d "$OUTPUT_DIR" ]; then
    echo "[WARN] Output directory $OUTPUT_DIR already exists. Remove it to re-run."
    echo "       To re-run: rm -rf $OUTPUT_DIR && bash picrust2_install_run.sh"
else
    echo "[INFO] Running PICRUSt2 pipeline (this takes ~60-90 min) ..."
    # --chunk_size 500: EPA-ng v0.3.8 silently fails with the default 5000
    # when the reference MSA is large (45 MB); 500 is safe and ~40 min for 7631 ASVs
    picrust2_pipeline.py \
        --study_fasta  "$INPUT_DIR/dna-sequences.fasta" \
        --input        "$INPUT_DIR/feature-table.biom" \
        --output       "$OUTPUT_DIR" \
        --stratified \
        --chunk_size   500 \
        --processes    "$THREADS" \
        2>&1 | tee picrust2_run.log
    echo "[INFO] Pipeline complete. Outputs in $OUTPUT_DIR/"
fi

# ---------------------------------------------------------------------------
# STEP 4 — Decompress key output files for R
# ---------------------------------------------------------------------------
echo "[INFO] Decompressing output files ..."

for f in \
    "$OUTPUT_DIR/pathways_out/path_abun_unstrat.tsv.gz" \
    "$OUTPUT_DIR/pathways_out/path_abun_strat.tsv.gz" \
    "$OUTPUT_DIR/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
    "$OUTPUT_DIR/EC_metagenome_out/pred_metagenome_strat.tsv.gz" \
    "$OUTPUT_DIR/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
    "$OUTPUT_DIR/KO_metagenome_out/pred_metagenome_strat.tsv.gz" \
    "$OUTPUT_DIR/intermediate/EC_predicted.tsv.gz" \
    "$OUTPUT_DIR/intermediate/KO_predicted.tsv.gz"
do
    if [ -f "$f" ]; then
        gunzip -kf "$f"
        echo "  decompressed: $f"
    fi
done

# ---------------------------------------------------------------------------
# STEP 5 — Summary
# ---------------------------------------------------------------------------
echo ""
echo "=== PICRUSt2 run complete ==="
echo "Key files for R analysis:"
echo "  Pathways (unstratified): $OUTPUT_DIR/pathways_out/path_abun_unstrat.tsv"
echo "  Pathways (stratified)  : $OUTPUT_DIR/pathways_out/path_abun_strat.tsv"
echo "  EC numbers (unstrat)   : $OUTPUT_DIR/EC_metagenome_out/pred_metagenome_unstrat.tsv"
echo "  EC numbers (strat)     : $OUTPUT_DIR/EC_metagenome_out/pred_metagenome_strat.tsv"
echo "  KO terms (unstrat)     : $OUTPUT_DIR/KO_metagenome_out/pred_metagenome_unstrat.tsv"
echo "  ASV-level EC predicted : $OUTPUT_DIR/intermediate/EC_predicted.tsv"
echo ""
echo "Next step: Rscript picrust2_analysis.R"
