import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import csv

from samples_array import samples


def calculate_molecular_weights(sample_list):
    results = []
    for s in sample_list:
        seq = s["sequence"].upper()
        bp = s["length_bp"]

        # Average dsDNA MW (~660 g/mol per bp, accounting for counterions)
        mw_avg_da = bp * 660.0
        mw_avg_kda = mw_avg_da / 1000.0
        mw_avg_mda = mw_avg_da / 1_000_000.0

        # Exact duplex MW based on monophosphate nucleotides:
        # Complementary strand has A->T, T->A, C->G, G->C
        nA = seq.count("A")
        nT = seq.count("T")
        nC = seq.count("C")
        nG = seq.count("G")
        # In duplex: total A-T pairs = nA + nT, total G-C pairs = nC + nG = bp
        # Average weight of dA-dT pair: 617.4 Da, dG-dC pair: 618.4 Da
        mw_exact_da = (nA + nT) * 617.4 + (nC + nG) * 618.4 + 157.9
        mw_exact_kda = mw_exact_da / 1000.0

        results.append({
            "sample_id": s["sample_id"],
            "start": s["start"],
            "end": s["end"],
            "length_bp": bp,
            "gc_content_pct": s["gc_content_pct"],
            "mw_avg_kda": round(mw_avg_kda, 2),
            "mw_avg_mda": round(mw_avg_mda, 3),
            "mw_exact_kda": round(mw_exact_kda, 2),
            "sequence": seq,
        })
    return results


def migration_distance_cm(bp, gel_length_cm=10.0, bp_min=100, bp_max=2000):

    log_bp = math.log10(max(50, bp))
    log_min = math.log10(bp_min)
    log_max = math.log10(bp_max)

    d_min = 1.5  # cm for 2000 bp
    d_max = 8.5  # cm for 100 bp

    # Linear semi-log relationship: d = d_max - slope * (log_bp - log_min)
    slope = (d_max - d_min) / (log_max - log_min)
    distance = d_max - slope * (log_bp - log_min)
    rf = distance / gel_length_cm
    return round(distance, 2), round(rf, 3)


def plot_simulation(fragment_data, output_file="gel_migration_simulation.png"):
    """
    Creates a dual-panel figure:
    Panel 1: Realistic UV Transilluminator Agarose Gel (Ladder + Pooled + Individual Lanes)
    Panel 2: Electrophoretic Mobility Semi-Log Calibration Curve (Distance vs log10(MW))
    """
    # Standard DNA Ladder (NEB/Thermo 100 bp - 1.5 kb reference)
    ladder_sizes = [1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100]

    plt.style.use("dark_background")
    fig = plt.figure(figsize=(15, 9), facecolor="#090d16")

    # Layout: GridSpec with 2 main subplots
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.25)
    ax_gel = fig.add_subplot(gs[0])
    ax_curve = fig.add_subplot(gs[1])

    # -------------------------------------------------------------
    # PANEL 1: SIMULATED AGAROSE GEL ELECTROPHORESIS
    # -------------------------------------------------------------
    ax_gel.set_facecolor("#060911")
    gel_height = 10.0  # cm
    num_lanes = 12
    lane_x = np.arange(1, num_lanes + 1)
    lane_names = ["Ladder", "Pool\n(All 10)"] + [f"F{i:02d}" for i in range(1, 11)]

    # Draw loading wells at top (y = 0.5 cm from top, represented at d = 0.5 cm)
    for lx in lane_x:
        well = patches.FancyBboxPatch(
            (lx - 0.35, 0.2), 0.7, 0.4,
            boxstyle="round,pad=0.04",
            facecolor="#182033", edgecolor="#334155", lw=1.2
        )
        ax_gel.add_patch(well)

    # Function to draw a realistic fluorescent DNA band
    def draw_band(ax, x_center, d_cm, bp, color="#f59e0b", width=0.68, height=0.14, alpha=0.9):
        # Draw multi-layered glow for realistic UV transilluminator fluorescence
        y_center = d_cm
        # Outer halo
        halo = patches.Rectangle(
            (x_center - width * 0.55, y_center - height * 1.2),
            width * 1.1, height * 2.4,
            facecolor=color, alpha=alpha * 0.18, edgecolor="none"
        )
        ax.add_patch(halo)
        # Inner glow
        mid = patches.Rectangle(
            (x_center - width * 0.45, y_center - height * 0.65),
            width * 0.9, height * 1.3,
            facecolor=color, alpha=alpha * 0.45, edgecolor="none"
        )
        ax.add_patch(mid)
        # Core intense band
        core = patches.Rectangle(
            (x_center - width * 0.38, y_center - height * 0.3),
            width * 0.76, height * 0.6,
            facecolor="#ffffff" if alpha > 0.8 else color, alpha=min(1.0, alpha * 0.9),
            edgecolor=color, lw=0.5
        )
        ax.add_patch(core)

    # 1. Draw Ladder bands (Lane 1)
    for sz in ladder_sizes:
        d, _ = migration_distance_cm(sz)
        is_landmark = sz in [500, 1000]
        alpha = 1.0 if is_landmark else 0.7
        draw_band(ax_gel, 1, d, sz, color="#38bdf8", alpha=alpha)
        # Label ladder bands on left margin
        ax_gel.text(
            0.42, d, f"{sz} bp",
            color="#7dd3fc" if not is_landmark else "#38bdf8",
            fontsize=8.5, va="center", ha="right",
            fontweight="bold" if is_landmark else "normal"
        )

    # 2. Draw Pooled Mixture (Lane 2)
    for frag in fragment_data:
        d, _ = migration_distance_cm(frag["length_bp"])
        draw_band(ax_gel, 2, d, frag["length_bp"], color="#34d399", alpha=0.85)

    # 3. Draw Individual Samples (Lanes 3..12)
    for i, frag in enumerate(fragment_data):
        lx = i + 3
        d, _ = migration_distance_cm(frag["length_bp"])
        draw_band(ax_gel, lx, d, frag["length_bp"], color="#fbbf24", alpha=0.92)
        # Display bp size above/below band
        ax_gel.text(
            lx, d - 0.22, f"{frag['length_bp']} bp",
            color="#fef08a", fontsize=7.2, ha="center", va="bottom", rotation=0
        )

    # Gel axes styling
    ax_gel.set_xlim(0.0, num_lanes + 0.9)
    ax_gel.set_ylim(gel_height + 0.2, 0.0)  # Invert so top is well (0 cm) and bottom is anode (10 cm)
    ax_gel.set_xticks(lane_x)
    ax_gel.set_xticklabels(lane_names, fontsize=8.5, fontweight="bold", color="#e2e8f0")
    ax_gel.set_ylabel("Migration Distance from Well (cm)", fontsize=11, fontweight="bold", color="#cbd5e1")
    ax_gel.set_title("A. Simulated Agarose Gel (1.2% TAE)\nEthidium Bromide / UV Transillumination",
                     fontsize=12, fontweight="bold", color="#f8fafc", pad=12)

    # Add Cathode (-) and Anode (+) Polarity indicators
    ax_gel.text(num_lanes + 0.5, 0.4, "(-) Cathode\n(Loading wells)", color="#ef4444", fontsize=8.5, ha="center", va="top", fontweight="bold")
    ax_gel.text(num_lanes + 0.5, gel_height - 0.3, "(+) Anode\n(Migration target)", color="#3b82f6", fontsize=8.5, ha="center", va="bottom", fontweight="bold")
    ax_gel.annotate("", xy=(num_lanes + 0.5, gel_height - 1.2), xytext=(num_lanes + 0.5, 1.2),
                    arrowprops=dict(arrowstyle="->", color="#94a3b8", lw=2.0))

    # Frame border for gel
    for spine in ax_gel.spines.values():
        spine.set_color("#334155")
        spine.set_linewidth(1.5)

    # -------------------------------------------------------------
    # PANEL 2: SEMI-LOG CALIBRATION CURVE (Distance vs log10(MW))
    # -------------------------------------------------------------
    ax_curve.set_facecolor("#0b1120")

    # Generate smooth theoretical standard line
    bp_line = np.linspace(100, 2000, 200)
    mw_line_kda = (bp_line * 660.0) / 1000.0
    d_line = [migration_distance_cm(b)[0] for b in bp_line]

    ax_curve.plot(mw_line_kda, d_line, color="#64748b", linestyle="--", lw=2, label="Theoretical Fit: $d = a - b \\cdot \\log_{10}(MW)$")

    # Plot ladder reference points
    ladder_mw = [(s * 660.0) / 1000.0 for s in ladder_sizes]
    ladder_d = [migration_distance_cm(s)[0] for s in ladder_sizes]
    ax_curve.scatter(ladder_mw, ladder_d, color="#38bdf8", s=50, zorder=4, label="DNA Ladder Markers", edgecolors="#ffffff", lw=0.8)

    # Plot the 10 DNA Fragments
    frag_mw = [f["mw_avg_kda"] for f in fragment_data]
    frag_d = [migration_distance_cm(f["length_bp"])[0] for f in fragment_data]

    scatter = ax_curve.scatter(
        frag_mw, frag_d,
        c=[f["length_bp"] for f in fragment_data],
        cmap="plasma_r",
        s=110,
        zorder=5,
        edgecolors="#ffffff",
        lw=1.2,
        label="10 Sample Fragments"
    )

    # Annotate each sample point
    for f, mw, d in zip(fragment_data, frag_mw, frag_d):
        ax_curve.annotate(
            f"{f['sample_id']} ({f['length_bp']} bp\n{mw:.1f} kDa)",
            xy=(mw, d),
            xytext=(mw * 1.08, d - 0.25 if mw < 500 else d + 0.25),
            fontsize=7.5,
            color="#f1f5f9",
            arrowprops=dict(arrowstyle="-|>", color="#94a3b8", lw=0.7, alpha=0.8),
            bbox=dict(boxstyle="round,pad=0.2", facecolor="#1e293b", edgecolor="#475569", alpha=0.85)
        )

    ax_curve.set_xscale("log")
    ax_curve.set_xlabel("Molecular Weight (kDa, log scale)", fontsize=11, fontweight="bold", color="#cbd5e1")
    ax_curve.set_ylabel("Migration Distance $d$ (cm)", fontsize=11, fontweight="bold", color="#cbd5e1")
    ax_curve.set_title("B. Electrophoretic Mobility Calibration Curve\nLogarithmic Molecular Weight vs. Migration Distance",
                       fontsize=12, fontweight="bold", color="#f8fafc", pad=12)
    ax_curve.grid(True, which="both", ls=":", color="#1e293b", alpha=0.7)

    # Secondary X-axis on top showing corresponding Base Pairs (bp)
    ax_top = ax_curve.secondary_xaxis(
        'top',
        functions=(lambda mw: (mw * 1000.0) / 660.0, lambda bp: (bp * 660.0) / 1000.0)
    )
    ax_top.set_xlabel("Double-Stranded DNA Size (bp)", fontsize=10, fontweight="bold", color="#93c5fd", labelpad=8)
    ax_top.tick_params(colors="#93c5fd")

    # Legend & annotations
    leg = ax_curve.legend(loc="lower left", facecolor="#111827", edgecolor="#374151", fontsize=9)
    for text in leg.get_texts():
        text.set_color("#f3f4f6")

    # Annotation box with physics explanation
    eq_text = (
        r"$\mathbf{Gel\ Physics:}$" "\n"
        r"$d \propto -\log_{10}(\mathrm{MW})$" "\n"
        r"Average dsDNA: $660\ \mathrm{Da/bp}$" "\n"
        r"Sieving: Smaller molecules migrate faster"
    )
    ax_curve.text(
        0.97, 0.97, eq_text,
        transform=ax_curve.transAxes,
        fontsize=8.5,
        va="top",
        ha="right",
        color="#e0e7ff",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#1e1b4b", edgecolor="#6366f1", alpha=0.9)
    )

    # Styling and padding adjustments
    for spine in ax_curve.spines.values():
        spine.set_color("#334155")
        spine.set_linewidth(1.5)

    fig.subplots_adjust(top=0.90, bottom=0.10, left=0.06, right=0.96, wspace=0.22)
    plt.savefig(output_file, dpi=250)
    plt.close()
    print(f"Generated simulation plot: {output_file}")


def export_results_csv(fragment_data, output_csv="gel_migration_results.csv"):
    """Export complete quantitative migration data."""
    fieldnames = [
        "sample_id",
        "length_bp",
        "mw_avg_kda",
        "mw_avg_mda",
        "mw_exact_kda",
        "gc_content_pct",
        "migration_distance_cm",
        "relative_mobility_rf",
    ]
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for f_data in fragment_data:
            d, rf = migration_distance_cm(f_data["length_bp"])
            row = {
                "sample_id": f_data["sample_id"],
                "length_bp": f_data["length_bp"],
                "mw_avg_kda": f_data["mw_avg_kda"],
                "mw_avg_mda": f_data["mw_avg_mda"],
                "mw_exact_kda": f_data["mw_exact_kda"],
                "gc_content_pct": f_data["gc_content_pct"],
                "migration_distance_cm": d,
                "relative_mobility_rf": rf,
            }
            writer.writerow(row)
    print(f"Exported migration metrics to: {output_csv}")


def main():
    print("Calculating molecular weights and electrophoretic mobilities...")
    frag_results = calculate_molecular_weights(samples)

    # Sort results by size descending to display like a real migration ladder
    frag_results_sorted = sorted(frag_results, key=lambda x: x["length_bp"], reverse=True)

    print("\n" + "=" * 92)
    print(f"{'Sample':<12} {'Length':<10} {'MW (kDa)':<12} {'MW (MDa)':<12} {'GC (%)':<9} {'Distance (cm)':<15} {'Mobility (Rf)':<12}")
    print("=" * 92)
    for f in frag_results_sorted:
        d, rf = migration_distance_cm(f["length_bp"])
        print(
            f"{f['sample_id']:<12} {f['length_bp']:<10} {f['mw_avg_kda']:<12.2f} "
            f"{f['mw_avg_mda']:<12.3f} {f['gc_content_pct']:<9.2f} {d:<15.2f} {rf:<12.3f}"
        )
    print("=" * 92)

    plot_simulation(frag_results, "gel_migration_simulation.png")
    export_results_csv(frag_results, "gel_migration_results.csv")


if __name__ == "__main__":
    main()
