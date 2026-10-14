#!/usr/bin/env bash
# ==============================================================================
# Complete Codon Analysis Pipeline Runner
# Location: /Users/toe/Documents/University/2026-2027/bioinfo/codon_analysis/
# ==============================================================================

set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$DIR"

echo "================================================================="
echo "🧬 BIOINFORMATICS CODON ANALYSIS PIPELINE"
echo "================================================================="
echo "Working Directory: $DIR"
echo ""

# 1. Paths
ECOLI_GENOME="/Users/toe/Documents/University/2026-2027/bioinfo/genome_datasets/ESCHERICHIA COLI/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
ECOLI_CDS="/Users/toe/Documents/University/2026-2027/bioinfo/genome_datasets/ESCHERICHIA COLI/ncbi_dataset/data/GCF_000005845.2/cds_from_genomic.fna"

echo "Step 1: Analyzing E. coli Whole Genome across ALL 6 frames..."
python3 codon_analyzer.py \
  --file "$ECOLI_GENOME" \
  --frame all \
  --top 10 \
  --chart "$DIR/ecoli_whole_genome_top10.png" \
  --csv "$DIR/ecoli_whole_genome_frequencies.csv"

echo ""
echo "Step 2: Analyzing E. coli Coding Sequences (CDS)..."
python3 codon_analyzer.py \
  --file "$ECOLI_CDS" \
  --frame cds \
  --top 10 \
  --chart "$DIR/ecoli_cds_top10.png" \
  --csv "$DIR/ecoli_cds_frequencies.csv"

echo ""
echo "Step 3: Generating LLM Prompts & Dietary Amino Acid Analysis..."
python3 llm_prompt_generator.py

echo ""
echo "================================================================="
echo "✅ PIPELINE COMPLETE!"
echo "================================================================="
echo "Generated Artifacts:"
echo "  📊 Charts:"
echo "     - $DIR/ecoli_whole_genome_top10.png"
echo "     - $DIR/ecoli_cds_top10.png"
echo "  💾 Data Tables:"
echo "     - $DIR/ecoli_whole_genome_frequencies.csv"
echo "     - $DIR/ecoli_cds_frequencies.csv"
echo "  🤖 LLM Prompts:"
echo "     - $DIR/prompt_whole_genome_glycine.txt"
echo "     - $DIR/prompt_cds_valine.txt"
echo ""
echo "To view charts on Mac:"
echo "  open ecoli_whole_genome_top10.png"
echo ""
echo "To launch the Interactive Web Dashboard:"
echo "  streamlit run app.py"
echo "================================================================="
