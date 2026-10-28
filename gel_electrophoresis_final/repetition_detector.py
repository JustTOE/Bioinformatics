#!/usr/bin/env python3
"""
Brute-Force DNA Repetition Detector & NCBI Sequence Analyzer
=============================================================
Detects repetitions between 6 and 10 nucleotides in DNA sequences,
produced by a brute-force engine generating nucleotide combinations
between 3 and 6 letters.

Supports:
- Local FASTA files (default: NZ_JBYXHI010000334.1.fasta)
- NCBI Entrez retrieval for arbitrary sequences (1000 - 3000 nt)
- Tandem repeat detection (SSRs / microsatellites)
- Interspersed repeated pattern detection (k-mers of length 6-10 with frequency >= 2)
- Visual distribution mapping and CSV/JSON export
"""

import argparse
import itertools
import json
import math
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from typing import Dict, List, Optional, Tuple


class BruteForceRepetitionEngine:
    """
    Bioinformatics engine for generating k-mer combinations and scanning DNA sequences
    for exact tandem repetitions and high-frequency repeated sequence patterns.
    """

    ALPHABET = ("A", "C", "G", "T")

    def __init__(self, sequence: str, header: str = "Custom Sequence"):
        self.raw_sequence = sequence.strip()
        self.sequence = "".join(self.raw_sequence.split()).upper()
        self.header = header
        self._validate_sequence()

    def _validate_sequence(self):
        """Ensure sequence contains valid IUPAC DNA nucleotides."""
        invalid_chars = set(self.sequence) - set(self.ALPHABET) - {"N"}
        if invalid_chars:
            raise ValueError(f"Sequence contains invalid characters: {invalid_chars}")

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def gc_content(self) -> float:
        if not self.sequence:
            return 0.0
        gc = self.sequence.count("G") + self.sequence.count("C")
        return round((gc / len(self.sequence)) * 100, 2)

    @classmethod
    def from_fasta_file(cls, filepath: str) -> "BruteForceRepetitionEngine":
        """Load sequence from a local FASTA file."""
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"FASTA file not found at: {filepath}")

        header = ""
        seq_parts = []
        with open(filepath, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if not header:
                        header = line[1:].strip()
                else:
                    seq_parts.append(line)

        return cls(sequence="".join(seq_parts), header=header)

    @classmethod
    def from_ncbi(cls, accession_id: str, email: str = "bioinfo_user@example.com") -> "BruteForceRepetitionEngine":
        """
        Fetch arbitrary nucleotide sequence from NCBI Entrez API.
        Validates whether sequence length falls within requested 1000 - 3000 nt range.
        """
        clean_id = accession_id.strip()
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "nuccore",
            "id": clean_id,
            "rettype": "fasta",
            "retmode": "text",
            "email": email,
        }
        query_string = urllib.parse.urlencode(params)
        req_url = f"{base_url}?{query_string}"

        req = urllib.request.Request(req_url, headers={"User-Agent": "NCBI-Sequence-Fetcher/1.0"})
        try:
            with urllib.request.urlopen(req, timeout=12) as response:
                content = response.read().decode("utf-8").strip()
        except urllib.error.URLError as e:
            raise ConnectionError(f"Failed to connect to NCBI Entrez: {e}")

        lines = content.splitlines()
        if not lines or not lines[0].startswith(">"):
            raise ValueError(f"Invalid response returned from NCBI for ID '{clean_id}':\n{content[:200]}")

        header = lines[0][1:].strip()
        sequence = "".join(line.strip() for line in lines[1:])

        print(f"Retrieved sequence from NCBI: {header}")
        print(f"Sequence length: {len(sequence)} nt")

        if len(sequence) < 1000 or len(sequence) > 3000:
            print(
                f"[Notice] Retrieved sequence length ({len(sequence)} nt) is outside the typical "
                f"1000-3000 nt window, but can still be processed."
            )

        return cls(sequence=sequence, header=header)

    @staticmethod
    def generate_combinations(min_k: int = 3, max_k: int = 6) -> Dict[int, List[str]]:
        """
        Brute-force generator of all nucleotide combinations for lengths min_k to max_k.
        For DNA (4 letters):
        k=3 -> 64 combinations
        k=4 -> 256 combinations
        k=5 -> 1024 combinations
        k=6 -> 4096 combinations
        Total = 5440 combinations.
        """
        combinations_by_k = {}
        for k in range(min_k, max_k + 1):
            combos = ["".join(p) for p in itertools.product(BruteForceRepetitionEngine.ALPHABET, repeat=k)]
            combinations_by_k[k] = combos
        return combinations_by_k

    def detect_tandem_repeats(
        self,
        min_motif: int = 3,
        max_motif: int = 6,
        min_total_len: int = 6,
        max_total_len: int = 10,
    ) -> List[Dict]:
        """
        Detects tandem repetitions where a combination of 3 to 6 letters
        is repeated consecutively r >= 2 times such that the total length
        is between min_total_len and max_total_len (6 to 10 letters).
        """
        combos_dict = self.generate_combinations(min_motif, max_motif)
        results = []
        seen_regions = set()

        # Iterate over each motif length k
        for k in range(min_motif, max_motif + 1):
            combos = combos_dict[k]
            # Determine valid repetition counts r such that min_total_len <= k * r <= max_total_len
            min_r = max(2, math.ceil(min_total_len / k))
            max_r = max_total_len // k

            if min_r > max_r:
                continue

            for motif in combos:
                for r in range(min_r, max_r + 1):
                    rep_pattern = motif * r
                    pattern_len = len(rep_pattern)

                    # Search for rep_pattern across sequence
                    start_idx = 0
                    while True:
                        idx = self.sequence.find(rep_pattern, start_idx)
                        if idx == -1:
                            break

                        start_1 = idx + 1
                        end_1 = idx + pattern_len
                        region_key = (start_1, end_1, rep_pattern)

                        if region_key not in seen_regions:
                            seen_regions.add(region_key)
                            gc = (rep_pattern.count("G") + rep_pattern.count("C")) / pattern_len * 100
                            results.append({
                                "type": "Tandem Repeat",
                                "motif": motif,
                                "motif_length": k,
                                "repeat_count": r,
                                "total_length": pattern_len,
                                "pattern": rep_pattern,
                                "start_position": start_1,
                                "end_position": end_1,
                                "gc_content_pct": round(gc, 2),
                            })
                        start_idx = idx + 1

        # Sort chronologically by start position, then by length
        results.sort(key=lambda x: (x["start_position"], -x["total_length"]))
        return results

    def detect_interspersed_repeats(
        self,
        min_len: int = 6,
        max_len: int = 10,
        min_occurrences: int = 2,
    ) -> List[Dict]:
        """
        Detects all repeated sequence patterns of length 6 to 10 letters
        that occur >= min_occurrences times anywhere across the sequence.
        """
        results = []
        for L in range(min_len, max_len + 1):
            counts = defaultdict(list)
            for i in range(len(self.sequence) - L + 1):
                sub = self.sequence[i : i + L]
                counts[sub].append(i + 1)  # 1-based start

            for pattern, positions in counts.items():
                if len(positions) >= min_occurrences:
                    gc = (pattern.count("G") + pattern.count("C")) / L * 100
                    results.append({
                        "type": "Interspersed Repeat",
                        "pattern": pattern,
                        "length": L,
                        "occurrences_count": len(positions),
                        "positions": positions,
                        "gc_content_pct": round(gc, 2),
                    })

        # Sort by occurrence count descending, then length descending
        results.sort(key=lambda x: (-x["occurrences_count"], -x["length"], x["pattern"]))
        return results

    def generate_visualization(self, tandem_repeats: List[Dict], output_image: str = "repetition_map.png") -> None:
        """
        Visualizes the spatial distribution and frequency of detected repetitions
        across the DNA sequence.
        """
        import matplotlib.pyplot as plt

        plt.style.use("dark_background")
        fig, (ax_map, ax_hist) = plt.subplots(
            2, 1, figsize=(13, 8), facecolor="#0b0f19", gridspec_kw={"height_ratios": [1.4, 1.0]}
        )

        # -------------------------------------------------------------
        # Subplot 1: Sequence Map (Lollipop / Density of Repetitions)
        # -------------------------------------------------------------
        ax_map.set_facecolor("#0f172a")
        seq_len = self.length

        # Draw sequence backbone
        ax_map.plot([1, seq_len], [0, 0], color="#475569", lw=4, zorder=1)

        # Color mapping by motif length
        color_map = {3: "#38bdf8", 4: "#34d399", 5: "#fbbf24", 6: "#f43f5e"}

        positions = [t["start_position"] for t in tandem_repeats]
        lengths = [t["total_length"] for t in tandem_repeats]
        motif_lens = [t["motif_length"] for t in tandem_repeats]

        for t in tandem_repeats:
            x = t["start_position"]
            y = t["total_length"]
            c = color_map.get(t["motif_length"], "#a855f7")
            ax_map.vlines(x, 0, y, color=c, alpha=0.6, lw=1.2, zorder=2)
            ax_map.scatter(x, y, color=c, s=36, zorder=3, edgecolors="#ffffff", lw=0.5)

        # Annotate longest repetitions (length 8, 9, 10)
        for t in tandem_repeats:
            if t["total_length"] >= 8:
                ax_map.annotate(
                    f"{t['pattern']}\n({t['total_length']}nt, pos {t['start_position']})",
                    xy=(t["start_position"], t["total_length"]),
                    xytext=(t["start_position"], t["total_length"] + 0.6),
                    fontsize=7.2,
                    ha="center",
                    color="#f8fafc",
                    arrowprops=dict(arrowstyle="->", color="#94a3b8", lw=0.8),
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="#1e293b", edgecolor="#475569", alpha=0.9),
                )

        ax_map.set_xlim(-20, seq_len + 50)
        ax_map.set_ylim(-0.5, 12.5)
        ax_map.set_xlabel("Nucleotide Position along Sequence (1 to N)", fontsize=10, fontweight="bold", color="#cbd5e1")
        ax_map.set_ylabel("Repetition Length (nt)", fontsize=10, fontweight="bold", color="#cbd5e1")
        ax_map.set_title(
            f"A. Genomic Distribution of Repetitions across Sequence\n({self.header[:70]}... | {seq_len} bp, {self.gc_content}% GC)",
            fontsize=11, fontweight="bold", color="#f8fafc", pad=10
        )
        ax_map.grid(True, ls=":", color="#334155", alpha=0.5)

        # Custom Legend for Motif Size
        for k in sorted(color_map.keys()):
            ax_map.scatter([], [], color=color_map[k], label=f"{k}-mer motif unit (k={k})", s=40)
        leg = ax_map.legend(loc="upper right", facecolor="#1e293b", edgecolor="#475569", fontsize=8.5)
        for text in leg.get_texts():
            text.set_color("#f1f5f9")

        # -------------------------------------------------------------
        # Subplot 2: Repetition Length Distribution Bar Chart
        # -------------------------------------------------------------
        ax_hist.set_facecolor("#0f172a")
        length_counts = defaultdict(int)
        for t in tandem_repeats:
            length_counts[t["total_length"]] += 1

        all_lengths = list(range(6, 11))
        counts = [length_counts[l] for l in all_lengths]

        bars = ax_hist.bar(all_lengths, counts, color="#6366f1", width=0.55, edgecolor="#a5b4fc", lw=1.2)
        for bar, count in zip(bars, counts):
            if count > 0:
                ax_hist.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 1.0,
                    f"{count}",
                    ha="center", va="bottom", fontsize=9, fontweight="bold", color="#e0e7ff"
                )

        ax_hist.set_xticks(all_lengths)
        ax_hist.set_xticklabels([f"{l} letters" for l in all_lengths], fontsize=9.5, color="#cbd5e1")
        ax_hist.set_xlabel("Repetition Length Window (6 to 10 letters)", fontsize=10, fontweight="bold", color="#cbd5e1")
        ax_hist.set_ylabel("Detected Occurrences", fontsize=10, fontweight="bold", color="#cbd5e1")
        ax_hist.set_ylim(0, max(counts) + 12 if counts else 10)
        ax_hist.set_title("B. Frequency of Repetitions by Total Length", fontsize=11, fontweight="bold", color="#f8fafc", pad=10)
        ax_hist.grid(axis="y", ls=":", color="#334155", alpha=0.5)

        plt.tight_layout()
        plt.savefig(output_image, dpi=220, bbox_inches="tight")
        plt.close()
        print(f"Generated visualization: {output_image}")


def export_csv(records: List[Dict], filepath: str) -> None:
    """Write records to CSV."""
    import csv
    if not records:
        return
    fieldnames = list(records[0].keys())
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    print(f"Exported {len(records)} records to CSV: {filepath}")


def export_json(records: List[Dict], filepath: str) -> None:
    """Write records to JSON."""
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(records, f, indent=2)
    print(f"Exported {len(records)} records to JSON: {filepath}")


def main():
    parser = argparse.ArgumentParser(
        description="Brute-Force DNA Repetition Detector (6-10 letters, 3-6 letter combinations)."
    )
    parser.add_argument(
        "--fasta",
        default="NZ_JBYXHI010000334.1.fasta",
        help="Path to local FASTA sequence file (default: NZ_JBYXHI010000334.1.fasta).",
    )
    parser.add_argument(
        "--ncbi-id",
        default=None,
        help="NCBI accession ID to fetch directly (e.g., NZ_JBYXHI010000334.1, NC_001416.1).",
    )
    parser.add_argument(
        "--min-motif",
        type=int,
        default=3,
        help="Minimum motif unit length in letters for brute force engine (default: 3).",
    )
    parser.add_argument(
        "--max-motif",
        type=int,
        default=6,
        help="Maximum motif unit length in letters for brute force engine (default: 6).",
    )
    parser.add_argument(
        "--min-len",
        type=int,
        default=6,
        help="Minimum total repetition length (default: 6).",
    )
    parser.add_argument(
        "--max-len",
        type=int,
        default=10,
        help="Maximum total repetition length (default: 10).",
    )
    parser.add_argument(
        "--mode",
        choices=["tandem", "interspersed", "both"],
        default="both",
        help="Repetition detection mode (default: both).",
    )
    parser.add_argument(
        "--output-csv",
        default="detected_repetitions.csv",
        help="Output CSV filepath.",
    )
    parser.add_argument(
        "--output-json",
        default="detected_repetitions.json",
        help="Output JSON filepath.",
    )
    parser.add_argument(
        "--plot",
        default="repetition_map.png",
        help="Output chart image path.",
    )

    args = parser.parse_args()

    # Load sequence
    if args.ncbi_id:
        print(f"Fetching accession '{args.ncbi_id}' from NCBI...")
        engine = BruteForceRepetitionEngine.from_ncbi(args.ncbi_id)
    else:
        print(f"Loading local sequence from '{args.fasta}'...")
        engine = BruteForceRepetitionEngine.from_fasta_file(args.fasta)

    print("\n" + "=" * 75)
    print(f"Sequence Header : {engine.header}")
    print(f"Sequence Length : {engine.length} nucleotides")
    print(f"GC Content      : {engine.gc_content}%")
    print("=" * 75)

    # Brute Force Engine Stats
    combos = engine.generate_combinations(args.min_motif, args.max_motif)
    total_combos = sum(len(c) for c in combos.values())
    print(f"\nBrute Force Engine generated {total_combos} motif combinations:")
    for k, v in combos.items():
        print(f"  - Length {k} letters: {len(v)} combinations (4^{k})")

    tandem_results = []
    interspersed_results = []

    if args.mode in ["tandem", "both"]:
        print(f"\nScanning for Tandem Repetitions ({args.min_len}-{args.max_len} letters)...")
        tandem_results = engine.detect_tandem_repeats(
            min_motif=args.min_motif,
            max_motif=args.max_motif,
            min_total_len=args.min_len,
            max_total_len=args.max_len,
        )
        print(f"-> Found {len(tandem_results)} tandem repetition occurrences.")

    if args.mode in ["interspersed", "both"]:
        print(f"\nScanning for High-Frequency Repeated Patterns ({args.min_len}-{args.max_len} letters)...")
        interspersed_results = engine.detect_interspersed_repeats(
            min_len=args.min_len,
            max_len=args.max_len,
            min_occurrences=2,
        )
        print(f"-> Found {len(interspersed_results)} unique repeated patterns appearing >= 2 times.")

    # Display Top Tandem Repeats
    if tandem_results:
        print("\n" + "-" * 75)
        print(f"{'Type':<15} {'Motif':<8} {'Repeat':<8} {'Length':<8} {'Positions':<14} {'Pattern':<15}")
        print("-" * 75)
        # Show repeats of length 8-10 first, then sample of length 6
        notable = [t for t in tandem_results if t["total_length"] >= 8]
        sample_6 = [t for t in tandem_results if t["total_length"] == 6][:10]

        for t in notable:
            print(
                f"{t['type']:<15} {t['motif']:<8} x{t['repeat_count']:<7} "
                f"{t['total_length']:<8} {t['start_position']}-{t['end_position']:<10} {t['pattern']:<15}"
            )
        if sample_6:
            print("... (sample of 6-letter tandem repeats):")
            for t in sample_6:
                print(
                    f"{t['type']:<15} {t['motif']:<8} x{t['repeat_count']:<7} "
                    f"{t['total_length']:<8} {t['start_position']}-{t['end_position']:<10} {t['pattern']:<15}"
                )
        print("-" * 75)

    # Display Top Interspersed Repeats
    if interspersed_results:
        print("\nTop 10 Most Frequent Repeated Sequences (Length 6-10 letters):")
        print("-" * 75)
        print(f"{'Pattern':<14} {'Length':<8} {'Occurrences':<13} {'GC (%)':<8} {'Sample Positions':<25}")
        print("-" * 75)
        for r in interspersed_results[:10]:
            pos_str = ", ".join(str(p) for p in r["positions"][:4])
            if len(r["positions"]) > 4:
                pos_str += "..."
            print(
                f"{r['pattern']:<14} {r['length']:<8} {r['occurrences_count']:<13} "
                f"{r['gc_content_pct']:<8.1f} {pos_str:<25}"
            )
        print("-" * 75)

    # Save to disk
    if tandem_results:
        export_csv(tandem_results, args.output_csv)
        export_json(tandem_results, args.output_json)
        engine.generate_visualization(tandem_results, args.plot)


if __name__ == "__main__":
    main()
