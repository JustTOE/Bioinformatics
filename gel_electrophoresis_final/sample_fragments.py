#!/usr/bin/env python3
"""
Fragment Sampler & Gel Electrophoresis Simulation
Samples random fragments from a DNA sequence within a specified length range
and exports both FASTA and CSV formats, with optional gel visualization.
"""

import argparse
import csv
import math
import os
import random
from typing import Dict, List, Tuple


def read_fasta(file_path: str) -> Tuple[str, str]:
    """Read first sequence and header from a FASTA file."""
    header = ""
    seq_parts = []
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if not header:
                    header = line[1:].strip()
            else:
                seq_parts.append(line)
    return header, "".join(seq_parts).upper()


def sample_fragments(
    sequence: str,
    num_samples: int = 10,
    min_len: int = 100,
    max_len: int = 3000,
    seed: int = 42,
) -> List[Dict]:
    """Sample random fragments of random lengths from the sequence."""
    if seed is not None:
        random.seed(seed)

    seq_len = len(sequence)
    effective_max_len = min(max_len, seq_len)

    if min_len > effective_max_len:
        raise ValueError(
            f"min_len ({min_len}) cannot exceed sequence length ({seq_len}) / max_len ({max_len})."
        )

    fragments = []
    for i in range(1, num_samples + 1):
        frag_len = random.randint(min_len, effective_max_len)
        start_0 = random.randint(0, seq_len - frag_len)
        end_0 = start_0 + frag_len
        subseq = sequence[start_0:end_0]

        # 1-based indexing for standard biological reporting
        start_1 = start_0 + 1
        end_1 = end_0

        gc_count = subseq.count("G") + subseq.count("C")
        gc_content = (gc_count / frag_len) * 100 if frag_len > 0 else 0.0

        fragments.append({
            "sample_id": f"Fragment_{i:02d}",
            "start": start_1,
            "end": end_1,
            "length_bp": frag_len,
            "gc_content_pct": round(gc_content, 2),
            "sequence": subseq,
        })

    return fragments


def export_fasta(fragments: List[Dict], output_path: str) -> None:
    """Export fragments to FASTA format."""
    with open(output_path, "w", encoding="utf-8") as f:
        for frag in fragments:
            header = (
                f">{frag['sample_id']} range={frag['start']}-{frag['end']} "
                f"length={frag['length_bp']}bp GC={frag['gc_content_pct']}%"
            )
            f.write(f"{header}\n")
            # Wrap at 70 characters per line
            seq = frag["sequence"]
            for j in range(0, len(seq), 70):
                f.write(f"{seq[j:j+70]}\n")


def export_csv(fragments: List[Dict], output_path: str) -> None:
    """Export fragment metadata and sequence to CSV."""
    fieldnames = ["sample_id", "start", "end", "length_bp", "gc_content_pct", "sequence"]
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for frag in fragments:
            writer.writerow(frag)


def simulate_gel(fragments: List[Dict], output_image: str) -> None:
    """
    Generate a simulated agarose gel electrophoresis visualization.
    Migration distance d is approximately inversely proportional to log10(bp).
    """
    import matplotlib.pyplot as plt
    import numpy as np

    # Standard DNA Ladder (e.g., 100 bp to 1500 bp ladder)
    ladder_sizes = [1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100]

    # Migration function: distance from well (top) to bottom
    # min_bp = 100 -> d ~ 0.9, max_bp = 2000 -> d ~ 0.15
    def calc_migration(bp: float) -> float:
        # Logarithmic migration model
        log_bp = math.log10(max(50, bp))
        log_max = math.log10(2000)
        log_min = math.log10(100)
        # Higher bp migrates less (top), lower bp migrates more (bottom)
        norm = (log_max - log_bp) / (log_max - log_min)
        return 0.15 + 0.75 * norm

    fig, ax = plt.subplots(figsize=(10, 8), facecolor="#1e1e1e")
    ax.set_facecolor("#121212")

    # Lanes: Lane 1 = Ladder, Lane 2 = Pooled Samples, Lanes 3..12 = Individual fragments
    num_lanes = 12
    lane_positions = list(range(1, num_lanes + 1))
    lane_labels = ["Ladder", "All (Pool)"] + [f"F{i}" for i in range(1, 11)]

    # Draw gel boundaries and wells
    for lane in lane_positions:
        # Well at the top
        well_rect = plt.Rectangle((lane - 0.35, 0.94), 0.7, 0.03, color="#2d3748", ec="#4a5568")
        ax.add_patch(well_rect)

    # Draw Ladder bands
    for size in ladder_sizes:
        y = 1.0 - calc_migration(size)
        # Intensity: standard reference bands like 500 bp or 1000 bp are brighter
        alpha = 0.95 if size in [500, 1000] else 0.75
        band = plt.Rectangle((1 - 0.3, y - 0.006), 0.6, 0.012, color="#60a5fa", alpha=alpha)
        ax.add_patch(band)
        ax.text(0.45, y, f"{size} bp", color="#93c5fd", fontsize=8, va="center", ha="right")

    # Draw Pooled Lane (Lane 2)
    for frag in fragments:
        y = 1.0 - calc_migration(frag["length_bp"])
        band = plt.Rectangle((2 - 0.3, y - 0.005), 0.6, 0.010, color="#34d399", alpha=0.7)
        ax.add_patch(band)

    # Draw Individual Lanes (Lanes 3..12)
    for idx, frag in enumerate(fragments):
        lane_idx = idx + 3
        y = 1.0 - calc_migration(frag["length_bp"])
        band = plt.Rectangle((lane_idx - 0.3, y - 0.006), 0.6, 0.012, color="#a78bfa", alpha=0.9)
        ax.add_patch(band)
        ax.text(lane_idx, y + 0.02, f"{frag['length_bp']} bp", color="#ddd6fe", fontsize=7.5, ha="center")

    ax.set_xlim(0, num_lanes + 0.8)
    ax.set_ylim(0.05, 1.05)
    ax.set_xticks(lane_positions)
    ax.set_xticklabels(lane_labels, color="#e2e8f0", fontsize=9, fontweight="bold")
    ax.set_yticks([])

    # Cathode (-) at top, Anode (+) at bottom
    ax.text(num_lanes + 0.5, 0.95, "(-) Cathode", color="#f87171", fontsize=9, ha="left", va="center", fontweight="bold")
    ax.text(num_lanes + 0.5, 0.10, "(+) Anode", color="#60a5fa", fontsize=9, ha="left", va="center", fontweight="bold")
    ax.annotate("", xy=(num_lanes + 0.5, 0.18), xytext=(num_lanes + 0.5, 0.88),
                arrowprops=dict(arrowstyle="->", color="#94a3b8", lw=1.5))

    ax.set_title(
        "Simulated Agarose Gel Electrophoresis (10 Sampled Fragments)\n"
        r"Migration $\propto -\log_{10}(\mathrm{bp})$",
        color="#f8fafc",
        fontsize=12,
        fontweight="bold",
        pad=15,
    )

    # Gel frame
    for spine in ax.spines.values():
        spine.set_color("#475569")
        spine.set_linewidth(1.5)

    plt.tight_layout()
    plt.savefig(output_image, dpi=200, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Sample random DNA fragments and simulate gel electrophoresis."
    )
    parser.add_argument(
        "--fasta",
        default="NZ_JBYXHI010000334.1.fasta",
        help="Path to input FASTA file.",
    )
    parser.add_argument(
        "--num-samples",
        type=int,
        default=10,
        help="Number of random samples to extract (default: 10).",
    )
    parser.add_argument(
        "--min-len",
        type=int,
        default=100,
        help="Minimum fragment length in bases (default: 100).",
    )
    parser.add_argument(
        "--max-len",
        type=int,
        default=3000,
        help="Maximum fragment length in bases (default: 3000).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42).",
    )
    parser.add_argument(
        "--output-fasta",
        default="sampled_fragments.fasta",
        help="Output FASTA file path.",
    )
    parser.add_argument(
        "--output-csv",
        default="sampled_fragments.csv",
        help="Output CSV file path.",
    )
    parser.add_argument(
        "--plot-gel",
        action="store_true",
        default=True,
        help="Generate simulated gel electrophoresis image (default: True).",
    )
    parser.add_argument(
        "--gel-image",
        default="gel_simulation.png",
        help="Output image file name for gel plot.",
    )

    args = parser.parse_args()

    if not os.path.exists(args.fasta):
        raise FileNotFoundError(f"FASTA file not found: {args.fasta}")

    header, seq = read_fasta(args.fasta)
    print(f"Loaded sequence: {header}")
    print(f"Sequence length: {len(seq)} bp")

    effective_max = min(args.max_len, len(seq))
    if effective_max < args.max_len:
        print(
            f"Note: Requested max length ({args.max_len} bp) exceeds sequence length "
            f"({len(seq)} bp). Sampling range adjusted to [{args.min_len}, {effective_max}] bp."
        )

    fragments = sample_fragments(
        sequence=seq,
        num_samples=args.num_samples,
        min_len=args.min_len,
        max_len=args.max_len,
        seed=args.seed,
    )

    export_fasta(fragments, args.output_fasta)
    export_csv(fragments, args.output_csv)
    print(f"Saved {len(fragments)} fragments to FASTA: {args.output_fasta}")
    print(f"Saved fragment details to CSV: {args.output_csv}")

    if args.plot_gel:
        simulate_gel(fragments, args.gel_image)
        print(f"Generated gel electrophoresis simulation plot: {args.gel_image}")

    print("\nSample Summary:")
    print("-" * 65)
    print(f"{'Sample ID':<13} {'Start':<8} {'End':<8} {'Length (bp)':<14} {'GC (%)':<8}")
    print("-" * 65)
    for frag in fragments:
        print(
            f"{frag['sample_id']:<13} {frag['start']:<8} {frag['end']:<8} "
            f"{frag['length_bp']:<14} {frag['gc_content_pct']:<8.2f}"
        )
    print("-" * 65)


if __name__ == "__main__":
    main()
