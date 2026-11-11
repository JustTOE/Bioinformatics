"""
visualizer.py - High-Resolution Genomic Visualization for Transposons & Intersections
======================================================================================
Provides publication-quality visual tracks, nucleotide sequence maps, and dot plots
comparing Ground Truth annotations with detection algorithms across intersected loci.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from typing import List, Dict, Any, Optional

# Color scheme for transposon families
TE_COLORS = {
    "TE1_Mariner": "#2563EB",   # Vibrant Blue
    "TE2_PiggyBac": "#DC2626",  # Red
    "TE3_Gypsy": "#16A34A",     # Green
    "TE4_Helitron": "#9333EA",  # Purple
    "Structural_TIR_Family": "#D97706", # Amber
    "Unknown": "#6B7280"        # Gray
}


def plot_genomic_tracks(simulation_result,
                        split_elements: Optional[List[Any]] = None,
                        naive_elements: Optional[List[Any]] = None,
                        denovo_elements: Optional[List[Any]] = None,
                        figsize: tuple = (14, 7)) -> plt.Figure:
    """
    Generate a publication-grade multi-track genomic visualization:
    - Track 1: Ground Truth (True insertion layout with fragments and nesting)
    - Track 2: Split-Aware Matcher Predictions (Resolved Compound Loci)
    - Track 3: Naive Continuous Matcher Predictions (Exposing Truncation/Failures)
    - Track 4: De Novo Structural Matcher (Optional, TIR/TSD discovery)
    """
    total_len = simulation_result.total_length

    # Determine number of tracks
    tracks = [("Ground Truth", "gt")]
    if split_elements is not None:
        tracks.append(("Split-Aware Matcher (Synteny Resolved)", "split"))
    if naive_elements is not None:
        tracks.append(("Naive Continuous Matcher (Truncated/Missed)", "naive"))
    if denovo_elements is not None:
        tracks.append(("Structural De Novo Matcher (TIR/TSD)", "denovo"))

    num_tracks = len(tracks)
    fig, ax = plt.subplots(figsize=figsize)

    ax.set_xlim(-5, total_len + 15)
    ax.set_ylim(-0.5, num_tracks * 1.5)

    # Invert Y so Track 1 is on top
    ax.invert_yaxis()

    track_height = 0.55

    for track_idx, (track_title, track_type) in enumerate(tracks):
        y_center = track_idx * 1.5 + 0.6

        # Background track bar (Host genomic DNA)
        ax.barh(y_center, total_len, left=0, height=0.10, color="#CBD5E1", edgecolor="#94A3B8", zorder=1)

        # Track label
        ax.text(-8, y_center, track_title, ha="right", va="center",
                fontsize=11, fontweight="bold", color="#1E293B")

        if track_type == "gt":
            # Plot Ground Truth Elements
            # To clearly visualize nesting, draw outer parents slightly lower or thicker
            gt_elements = sorted(simulation_result.elements.values(), key=lambda e: e.full_end - e.full_start, reverse=True)
            for elem in gt_elements:
                color = TE_COLORS.get(elem.family_name, "#2563EB")
                
                # Draw outer span bracket/dashed line if interrupted
                if elem.is_interrupted:
                    span_start = elem.full_start
                    span_width = elem.full_end - elem.full_start + 1
                    rect_dashed = patches.Rectangle(
                        (span_start, y_center - track_height / 2),
                        span_width, track_height,
                        fill=False, linestyle="--", edgecolor=color, linewidth=1.5, zorder=2
                    )
                    ax.add_patch(rect_dashed)

                # Draw individual physical fragments
                for frag in elem.fragments:
                    f_width = frag.end - frag.start + 1
                    rect = patches.Rectangle(
                        (frag.start, y_center - (track_height - 0.1) / 2),
                        f_width, track_height - 0.1,
                        facecolor=color, edgecolor="#0F172A", alpha=0.90, linewidth=1.2, zorder=3
                    )
                    ax.add_patch(rect)
                    # Label inside fragment if wide enough
                    if f_width >= 18:
                        lbl = f"{elem.element_id.split('_')[0]}\n({frag.fragment_label.split('_')[0]})"
                        ax.text(frag.start + f_width / 2, y_center, lbl,
                                ha="center", va="center", color="white", fontsize=8, fontweight="bold", zorder=4)

                # TSD markers
                if elem.tsd_start_5 is not None:
                    ax.plot(elem.tsd_start_5, y_center, 'o', color="#F59E0B", markersize=6, zorder=5)
                if elem.tsd_start_3 is not None:
                    ax.plot(elem.tsd_start_3, y_center, 'o', color="#F59E0B", markersize=6, zorder=5)

        elif track_type == "split":
            for elem in split_elements:
                color = TE_COLORS.get(elem.family_name, "#2563EB")
                
                # Outer span
                if elem.is_interrupted:
                    span_start = elem.start
                    span_width = elem.end - elem.start + 1
                    rect_dashed = patches.Rectangle(
                        (span_start, y_center - track_height / 2),
                        span_width, track_height,
                        fill=False, linestyle="-.", edgecolor=color, linewidth=1.8, zorder=2
                    )
                    ax.add_patch(rect_dashed)

                # Fragments
                for frag in elem.fragments:
                    f_width = frag.end - frag.start + 1
                    rect = patches.Rectangle(
                        (frag.start, y_center - (track_height - 0.1) / 2),
                        f_width, track_height - 0.1,
                        facecolor=color, edgecolor="#0F172A", alpha=0.85, linewidth=1.2, zorder=3
                    )
                    ax.add_patch(rect)
                    if f_width >= 18:
                        ax.text(frag.start + f_width / 2, y_center, f"{elem.element_id.split('_')[0]}",
                                ha="center", va="center", color="white", fontsize=8, fontweight="bold", zorder=4)

        elif track_type == "naive":
            for elem in naive_elements:
                color = TE_COLORS.get(elem.family_name, "#2563EB")
                e_width = elem.end - elem.start + 1
                rect = patches.Rectangle(
                    (elem.start, y_center - (track_height - 0.1) / 2),
                    e_width, track_height - 0.1,
                    facecolor=color, edgecolor="#DC2626" if elem.confidence_score < 0.8 else "#0F172A",
                    alpha=0.75, linewidth=1.5, zorder=3
                )
                ax.add_patch(rect)
                if e_width >= 18:
                    lbl = f"{elem.element_id.split('_')[0]}"
                    if elem.confidence_score < 0.8:
                        lbl += "\n(TRUNC)"
                    ax.text(elem.start + e_width / 2, y_center, lbl,
                            ha="center", va="center", color="white", fontsize=8, fontweight="bold", zorder=4)

        elif track_type == "denovo":
            for elem in denovo_elements:
                color = "#D97706"
                e_width = elem.end - elem.start + 1
                rect = patches.Rectangle(
                    (elem.start, y_center - (track_height - 0.1) / 2),
                    e_width, track_height - 0.1,
                    facecolor=color, edgecolor="#78350F", alpha=0.8, linewidth=1.2, zorder=3
                )
                ax.add_patch(rect)
                if e_width >= 20:
                    ax.text(elem.start + e_width / 2, y_center, "TIR Pair",
                            ha="center", va="center", color="white", fontsize=8, fontweight="bold", zorder=4)

    # Aesthetic styling
    ax.set_xlabel("Genomic Position (Base Pairs, 0-indexed)", fontsize=11, fontweight="bold", labelpad=10)
    ax.set_xticks(range(0, total_len + 25, 25))
    ax.tick_params(axis='x', which='both', labelsize=9)
    ax.yaxis.set_visible(False)
    ax.grid(axis='x', linestyle=':', alpha=0.4, color='#94A3B8')

    # Legend
    legend_handles = [
        patches.Patch(facecolor=TE_COLORS["TE1_Mariner"], label="TE1_Mariner (DNA/TIR)"),
        patches.Patch(facecolor=TE_COLORS["TE2_PiggyBac"], label="TE2_PiggyBac (DNA/TIR)"),
        patches.Patch(facecolor=TE_COLORS["TE3_Gypsy"], label="TE3_Gypsy (Retro/LTR)"),
        patches.Patch(facecolor=TE_COLORS["TE4_Helitron"], label="TE4_Helitron (DNA/TIR)"),
        patches.Patch(facecolor="#CBD5E1", edgecolor="#94A3B8", label="Host Genomic DNA"),
        patches.Patch(facecolor="none", edgecolor="#2563EB", linestyle="--", label="Interrupted Element Span")
    ]
    ax.legend(handles=legend_handles, loc="upper right", bbox_to_anchor=(1.0, 1.25),
              ncol=3, fontsize=9, frameon=True)

    plt.title(f"Transposable Element Architecture & Algorithmic Detection Map\n"
              f"Length: {total_len} bp | Scenario: {simulation_result.scenario_name}",
              fontsize=13, fontweight="bold", pad=28)

    plt.tight_layout()
    return fig


def plot_dotplot_synteny(query_seq: str, ref_library: Dict[str, str],
                         figsize: tuple = (10, 6), kmer: int = 5) -> plt.Figure:
    """
    Generate a 2D Dot-Plot showing query coordinate vs reference TE coordinates.
    Nested/interrupted transposons produce a characteristic stepped diagonal signature.
    """
    fig, axes = plt.subplots(1, len(ref_library), figsize=figsize, sharey=False)
    if len(ref_library) == 1:
        axes = [axes]

    query_len = len(query_seq)

    for ax, (te_name, ref_seq) in zip(axes, ref_library.items()):
        ref_len = len(ref_seq)
        dot_matrix = np.zeros((ref_len, query_len), dtype=bool)

        for i in range(ref_len - kmer + 1):
            sub_r = ref_seq[i:i + kmer]
            for j in range(query_len - kmer + 1):
                if query_seq[j:j + kmer] == sub_r:
                    dot_matrix[i, j] = True

        ax.imshow(dot_matrix, cmap="Blues", origin="lower", aspect="auto")
        ax.set_title(te_name.split("_")[0], fontsize=10, fontweight="bold")
        ax.set_xlabel("Query DNA (bp)", fontsize=9)
        ax.set_ylabel("Reference (bp)", fontsize=9)

    plt.suptitle("Synteny Dot-Plots (Stepped Diagonal Shows Interrupted/Nested Elements)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    return fig


def generate_html_sequence_viewer(dna_sequence: str,
                                  ground_truth_elements: Dict[str, Any],
                                  chars_per_line: int = 60) -> str:
    """
    Generate an interactive, beautifully highlighted HTML nucleotide sequence view.
    Each base is color-coded by its transposable element membership or host origin.
    """
    total_len = len(dna_sequence)
    # Map each position to metadata
    pos_meta = [None] * total_len

    for elem_id, elem in ground_truth_elements.items():
        color = TE_COLORS.get(elem.family_name, "#2563EB")
        for frag in elem.fragments:
            for p in range(frag.start, frag.end + 1):
                pos_meta[p] = {
                    "color": color,
                    "element": elem_id,
                    "fragment": frag.fragment_label,
                    "nested_level": elem.nesting_level
                }

    html = [
        """
        <div style="font-family: 'JetBrains Mono', 'Courier New', monospace; 
                    font-size: 13px; line-height: 1.8; background-color: #0F172A; 
                    padding: 20px; border-radius: 8px; color: #F8FAFC; overflow-x: auto;">
        <div style="margin-bottom: 15px; font-weight: bold; border-bottom: 1px solid #334155; padding-bottom: 8px;">
            <span>🧬 Nucleotide Sequence Inspection (200-500 bp)</span>
            <div style="float: right; font-size: 11px;">
                <span style="display:inline-block; width:12px; height:12px; background:#2563EB; margin-right:4px;"></span>TE1
                <span style="display:inline-block; width:12px; height:12px; background:#DC2626; margin-left:8px; margin-right:4px;"></span>TE2
                <span style="display:inline-block; width:12px; height:12px; background:#16A34A; margin-left:8px; margin-right:4px;"></span>TE3
                <span style="display:inline-block; width:12px; height:12px; background:#9333EA; margin-left:8px; margin-right:4px;"></span>TE4
                <span style="display:inline-block; width:12px; height:12px; background:#475569; margin-left:8px; margin-right:4px;"></span>Host
            </div>
        </div>
        """
    ]

    for line_start in range(0, total_len, chars_per_line):
        line_end = min(line_start + chars_per_line, total_len)
        chunk = dna_sequence[line_start:line_end]

        html.append(f'<div style="margin-bottom: 6px;">')
        html.append(f'<span style="color: #64748B; width: 50px; display: inline-block;">{line_start:04d}</span> ')

        for idx, base in enumerate(chunk):
            pos = line_start + idx
            m = pos_meta[pos]
            if m:
                bg = m["color"]
                title = f"Pos {pos}: {base} ({m['element']} - {m['fragment']}, Nesting: L{m['nested_level']})"
                html.append(
                    f'<span title="{title}" style="background-color: {bg}; color: #FFFFFF; '
                    f'padding: 1px 3px; border-radius: 2px; font-weight: bold; cursor: pointer;">{base}</span>'
                )
            else:
                title = f"Pos {pos}: {base} (Host Background DNA)"
                html.append(
                    f'<span title="{title}" style="color: #94A3B8; padding: 1px 3px;">{base}</span>'
                )

        html.append(f' <span style="color: #64748B; margin-left: 10px;">{line_end - 1:04d}</span>')
        html.append('</div>')

    html.append('</div>')
    return "".join(html)
