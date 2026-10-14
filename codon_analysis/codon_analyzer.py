"""
codon_analyzer.py - Codon Frequency Analysis and Visualization Engine
Calculates codon frequencies across all reading frames (+1, +2, +3, -1, -2, -3, or combined),
ranks codons from highest to lowest frequency, and creates a chart of the top 10 codons.
"""

import argparse
import os
import sys
from collections import Counter
from typing import Dict, List, Optional, Tuple, Union

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend suitable for scripts and servers
import matplotlib.pyplot as plt
import numpy as np

# Import translation functions
try:
    from genetic_code import (
        codon_to_amino_acid,
        get_codon_info,
        reverse_complement,
        AMINO_ACID_PROPERTIES,
    )
except ImportError:
    from .genetic_code import (
        codon_to_amino_acid,
        get_codon_info,
        reverse_complement,
        AMINO_ACID_PROPERTIES,
    )


# Color map for amino acid chemical classes
CLASS_COLORS = {
    'Hydrophobic / Nonpolar': '#2b5c8f',       # Blue
    'Hydrophobic / Start': '#1f77b4',          # Dark Blue
    'Hydrophobic / Aromatic': '#3470a3',       # Slate Blue
    'Polar / Uncharged': '#2ca02c',            # Green
    'Polar / Aromatic': '#3dbd3d',             # Light Green
    'Polar / Special (Disulfide)': '#17becf',  # Cyan
    'Special / Smallest': '#bcbd22',           # Olive
    'Special / Cyclic': '#e377c2',             # Pink
    'Acidic / Negatively Charged': '#d62728',  # Red
    'Basic / Positively Charged': '#ff7f0e',   # Orange
    'Basic / Aromatic': '#e6550d',             # Dark Orange
    'Stop Codon': '#7f7f7f',                   # Grey
    'Ambiguous': '#c7c7c7'                     # Light grey
}


def read_fasta(file_path: str) -> List[Tuple[str, str]]:
    """
    Parses a FASTA file and returns a list of (header, sequence) tuples.
    Supports both single-sequence genomes and multi-entry CDS files.
    """
    records = []
    current_header = None
    current_seq_chunks = []
    
    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(current_seq_chunks).upper()))
                    current_seq_chunks = []
                current_header = line[1:]
            else:
                current_seq_chunks.append(line)
        if current_header is not None:
            records.append((current_header, "".join(current_seq_chunks).upper()))
            
    return records


def extract_codons_from_sequence(seq: str, frame: int) -> List[str]:
    """
    Extracts all 3-nucleotide triplets from a sequence in a specified reading frame.
    
    Parameters:
        seq (str): Nucleotide sequence (DNA).
        frame (int): Frame: +1, +2, +3 for forward; -1, -2, -3 for reverse strand.
        
    Returns:
        List[str]: List of 3-nucleotide codons.
    """
    dna = seq.upper().replace("U", "T")
    
    if frame < 0:
        dna = reverse_complement(dna)
        offset = abs(frame) - 1
    else:
        offset = frame - 1
        
    codons = []
    seq_len = len(dna)
    for i in range(offset, seq_len - 2, 3):
        triplet = dna[i:i+3]
        # Keep only standard 3-base triplets
        if len(triplet) == 3 and all(b in "ATCG" for b in triplet):
            codons.append(triplet)
            
    return codons


def analyze_frequencies(
    fasta_path: str,
    frame_mode: str = "all",
    table: int = 1
) -> Dict[str, Union[int, float, List[Dict]]]:
    """
    Computes codon counts and frequencies across the requested reading frame(s).
    
    Parameters:
        fasta_path (str): Path to FASTA file.
        frame_mode (str): One of:
            'all'      -> All 6 reading frames (+1, +2, +3, -1, -2, -3) combined
            'forward'  -> Forward frames (+1, +2, +3) combined
            '+1', '+2', '+3', '-1', '-2', '-3' -> Individual specific frame
            'cds'      -> Frame +1 for all coding sequences (for CDS FASTA files)
        table (int): NCBI genetic code table (1 = Standard, 11 = Bacterial).
        
    Returns:
        Dict containing:
            - total_codons: int
            - total_nucleotides: int
            - gc_content_percent: float
            - frame_mode: str
            - ranked_codons: List[Dict] sorted from highest frequency to lowest
    """
    records = read_fasta(fasta_path)
    if not records:
        raise ValueError(f"No FASTA records found in {fasta_path}")
        
    total_nt = sum(len(seq) for _, seq in records)
    g_count = sum(seq.count("G") for _, seq in records)
    c_count = sum(seq.count("C") for _, seq in records)
    gc_percent = ((g_count + c_count) / total_nt * 100) if total_nt > 0 else 0.0
    
    codon_counter = Counter()
    
    if frame_mode.lower() == "all":
        frames_to_run = [+1, +2, +3, -1, -2, -3]
    elif frame_mode.lower() == "forward":
        frames_to_run = [+1, +2, +3]
    elif frame_mode.lower() == "reverse":
        frames_to_run = [-1, -2, -3]
    elif frame_mode in ["+1", "1"]:
        frames_to_run = [+1]
    elif frame_mode in ["+2", "2"]:
        frames_to_run = [+2]
    elif frame_mode in ["+3", "3"]:
        frames_to_run = [+3]
    elif frame_mode in ["-1"]:
        frames_to_run = [-1]
    elif frame_mode in ["-2"]:
        frames_to_run = [-2]
    elif frame_mode in ["-3"]:
        frames_to_run = [-3]
    elif frame_mode.lower() == "cds":
        # For CDS files, evaluate frame +1 for each CDS entry
        frames_to_run = [+1]
    else:
        frames_to_run = [+1, +2, +3, -1, -2, -3]
        
    for _, seq in records:
        for f in frames_to_run:
            codons = extract_codons_from_sequence(seq, f)
            codon_counter.update(codons)
            
    total_codons = sum(codon_counter.values())
    
    ranked_codons = []
    # Sort by count descending
    for rank, (codon, count) in enumerate(codon_counter.most_common(), start=1):
        freq_pct = (count / total_codons * 100) if total_codons > 0 else 0.0
        freq_per_thousand = (count / total_codons * 1000) if total_codons > 0 else 0.0
        info = get_codon_info(codon, table=table)
        
        ranked_codons.append({
            "rank": rank,
            "codon": codon,
            "codon_rna": codon.replace("T", "U"),
            "amino_acid": info["amino_acid_1"],
            "amino_acid_3": info["amino_acid_3"],
            "amino_acid_name": info["amino_acid_name"],
            "category": info["category"],
            "count": count,
            "frequency_percent": round(freq_pct, 4),
            "frequency_per_thousand": round(freq_per_thousand, 2),
            "is_start": info["is_start"],
            "is_stop": info["is_stop"]
        })
        
    return {
        "num_records": len(records),
        "total_nucleotides": total_nt,
        "gc_content_percent": round(gc_percent, 2),
        "total_codons": total_codons,
        "frame_mode": frame_mode,
        "ranked_codons": ranked_codons
    }


def plot_top_10_codons(
    ranked_codons: List[Dict],
    title: str = "Top 10 Most Frequent Codons",
    subtitle: Optional[str] = None,
    output_path: str = "top_10_codons.png"
) -> str:
    """
    Renders a publication-ready horizontal or vertical bar chart of the Top 10 most frequent codons.
    """
    top_10 = ranked_codons[:10]
    if not top_10:
        raise ValueError("Cannot plot: No codon data provided.")
        
    codons = [item["codon"] for item in top_10]
    labels = [f"{item['codon']}\n({item['amino_acid_3']} / {item['amino_acid']})" for item in top_10]
    frequencies = [item["frequency_percent"] for item in top_10]
    counts = [item["count"] for item in top_10]
    categories = [item["category"] for item in top_10]
    
    colors = [CLASS_COLORS.get(cat, '#2b5c8f') for cat in categories]
    
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, ax = plt.subplots(figsize=(11, 6.5), dpi=300)
    
    bars = ax.bar(labels, frequencies, color=colors, edgecolor='#333333', linewidth=1.1, width=0.65, zorder=3)
    
    # Add value annotations above bars
    max_freq = max(frequencies) if frequencies else 1.0
    ax.set_ylim(0, max_freq * 1.18)
    
    for bar, freq, count in zip(bars, frequencies, counts):
        height = bar.get_height()
        ax.annotate(
            f"{freq:.2f}%\n({count:,})",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 5),
            textcoords="offset points",
            ha='center', va='bottom',
            fontsize=9.5, fontweight='bold',
            color='#1a1a1a'
        )
        
    ax.set_title(title, fontsize=15, fontweight='bold', pad=18, color='#111111')
    if subtitle:
        ax.text(0.5, 1.02, subtitle, transform=ax.transAxes, ha='center', fontsize=10.5, color='#555555')
        
    ax.set_xlabel("Codon (Amino Acid 3-letter / 1-letter)", fontsize=11, fontweight='semibold', labelpad=10)
    ax.set_ylabel("Frequency (%) across reading frames", fontsize=11, fontweight='semibold', labelpad=10)
    
    # Custom Legend for Amino Acid Classes
    unique_cats = list(dict.fromkeys(categories))
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=CLASS_COLORS.get(cat, '#2b5c8f'), edgecolor='#333333')
        for cat in unique_cats
    ]
    ax.legend(
        legend_handles, unique_cats,
        title="Amino Acid Chemical Class",
        title_fontsize=9.5,
        fontsize=8.5,
        loc="upper right",
        frameon=True,
        facecolor="white",
        edgecolor="#cccccc"
    )
    
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return os.path.abspath(output_path)


def save_frequencies_to_csv(ranked_codons: List[Dict], csv_path: str):
    """
    Exports the complete sorted codon frequencies to a CSV file.
    """
    import csv
    os.makedirs(os.path.dirname(os.path.abspath(csv_path)), exist_ok=True)
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "rank", "codon", "codon_rna", "amino_acid", "amino_acid_3",
            "amino_acid_name", "category", "count", "frequency_percent",
            "frequency_per_thousand", "is_start", "is_stop"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(ranked_codons)


def main():
    parser = argparse.ArgumentParser(
        description="Codon Frequency Analyzer - Analyzes all reading frames and plots top 10 codons."
    )
    parser.add_argument(
        "--file", "-f",
        default="/Users/toe/Documents/University/2026-2027/bioinfo/genome_datasets/ESCHERICHIA COLI/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna",
        help="Path to FASTA genome or CDS file."
    )
    parser.add_argument(
        "--frame", "-r",
        default="all",
        choices=["all", "forward", "reverse", "+1", "+2", "+3", "-1", "-2", "-3", "cds"],
        help="Reading frame mode: 'all' (all 6 frames), 'forward' (+1,+2,+3), 'cds', or specific frame."
    )
    parser.add_argument(
        "--table", "-t",
        type=int,
        default=1,
        help="NCBI translation table (1 = Standard, 11 = Bacterial)."
    )
    parser.add_argument(
        "--top", "-n",
        type=int,
        default=10,
        help="Number of top codons to display in summary."
    )
    parser.add_argument(
        "--chart", "-c",
        default="top_10_codons.png",
        help="Output image path for the top 10 codons chart."
    )
    parser.add_argument(
        "--csv",
        default="codon_frequencies.csv",
        help="Output CSV path for full sorted codon table."
    )
    
    args = parser.parse_args()
    
    print(f"\n=================================================================")
    print(f"🧬 CODON FREQUENCY ANALYZER")
    print(f"=================================================================")
    print(f"Loading input file: {args.file}")
    
    if not os.path.exists(args.file):
        print(f"❌ Error: File not found at {args.file}")
        sys.exit(1)
        
    print(f"Analyzing reading frame(s): {args.frame.upper()} | Genetic code table: {args.table}...")
    res = analyze_frequencies(args.file, frame_mode=args.frame, table=args.table)
    
    print(f"\n📊 Genome / Dataset Overview:")
    print(f"  • Total FASTA records : {res['num_records']}")
    print(f"  • Total nucleotides   : {res['total_nucleotides']:,} bp")
    print(f"  • GC content          : {res['gc_content_percent']}%")
    print(f"  • Total codons counted: {res['total_codons']:,}")
    print(f"  • Reading frame mode  : {res['frame_mode']}")
    
    ranked = res["ranked_codons"]
    print(f"\n🏆 Top {args.top} Most Frequent Codons (Highest to Lowest):")
    print(f"{'Rank':<5} | {'Codon':<6} | {'AA':<4} | {'Name':<15} | {'Class':<22} | {'Count':<10} | {'Freq (%)':<9} | {'Per 1k':<7}")
    print("-" * 90)
    for row in ranked[:args.top]:
        print(
            f"{row['rank']:<5} | {row['codon']:<6} | {row['amino_acid_3']:<4} | {row['amino_acid_name']:<15} | "
            f"{row['category'][:22]:<22} | {row['count']:<10,d} | {row['frequency_percent']:<9.3f} | {row['frequency_per_thousand']:<7.1f}"
        )
        
    # Generate and save chart
    chart_title = f"Top {args.top} Most Frequent Codons ({res['frame_mode'].upper()} Frame(s))"
    sub_title = f"Total Codons: {res['total_codons']:,} | GC Content: {res['gc_content_percent']}%"
    chart_path = plot_top_10_codons(ranked, title=chart_title, subtitle=sub_title, output_path=args.chart)
    print(f"\n📈 Chart successfully generated: {chart_path}")
    
    # Save CSV
    save_frequencies_to_csv(ranked, args.csv)
    print(f"💾 Full 64-codon frequency table saved to: {os.path.abspath(args.csv)}")
    print(f"=================================================================\n")


if __name__ == "__main__":
    main()
