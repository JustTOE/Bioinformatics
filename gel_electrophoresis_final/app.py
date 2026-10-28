"""
app.py - Interactive Streamlit Application for DNA Repetition Detection
=======================================================================
Detects DNA repetitions between 6 and 10 nucleotides in DNA sequences,
produced by a brute-force combination engine (3 to 6 letters).
Supports local FASTA and direct NCBI Entrez retrieval.

Run with:
    streamlit run app.py
"""

import io
import json
import os
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

from repetition_detector import BruteForceRepetitionEngine

# Page configuration
st.set_page_config(
    page_title="DNA Repetition Detector",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom Styling
st.markdown(
    """
    <style>
    .metric-card {
        background-color: #1e293b;
        padding: 15px;
        border-radius: 8px;
        border: 1px solid #334155;
        text-align: center;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    .stTabs [data-baseweb="tab"] {
        background-color: #0f172a;
        border-radius: 6px;
        padding: 8px 16px;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# Header
st.title("🧬 DNA Repetition Detector & NCBI Analyzer")
st.markdown(
    "Detects **repetitions between 6 and 10 letters** produced by a **brute-force combination engine "
    "(3 to 6 letters)** in DNA sequences (1,000–3,000 nucleotides)."
)

# Sidebar: Controls & Input
st.sidebar.header("⚙️ Sequence Source")

input_option = st.sidebar.radio(
    "Select Input Method:",
    [
        "Local Sequence (NZ_JBYXHI010000334.1.fasta)",
        "Fetch from NCBI Entrez",
        "Upload FASTA File",
        "Paste DNA Sequence",
    ],
    index=0,
)

engine = None
header = ""
sequence = ""

if input_option == "Local Sequence (NZ_JBYXHI010000334.1.fasta)":
    local_path = "NZ_JBYXHI010000334.1.fasta"
    if os.path.exists(local_path):
        engine = BruteForceRepetitionEngine.from_fasta_file(local_path)
    else:
        st.error(f"Local file '{local_path}' not found.")

elif input_option == "Fetch from NCBI Entrez":
    ncbi_id = st.sidebar.text_input("NCBI Accession ID:", value="NZ_JBYXHI010000334.1")
    fetch_btn = st.sidebar.button("📥 Fetch Sequence")
    if fetch_btn or "ncbi_engine" in st.session_state:
        if fetch_btn:
            try:
                with st.spinner(f"Fetching '{ncbi_id}' from NCBI Entrez..."):
                    engine = BruteForceRepetitionEngine.from_ncbi(ncbi_id)
                    st.session_state["ncbi_engine"] = engine
            except Exception as e:
                st.sidebar.error(f"Error fetching from NCBI: {e}")
        else:
            engine = st.session_state.get("ncbi_engine")

elif input_option == "Upload FASTA File":
    uploaded = st.sidebar.file_uploader("Upload FASTA (.fasta, .fna, .txt)", type=["fasta", "fna", "txt"])
    if uploaded:
        content = uploaded.getvalue().decode("utf-8")
        lines = content.strip().splitlines()
        h = lines[0][1:].strip() if lines[0].startswith(">") else "Uploaded Sequence"
        s = "".join(l.strip() for l in lines if not l.startswith(">"))
        engine = BruteForceRepetitionEngine(s, header=h)

elif input_option == "Paste DNA Sequence":
    pasted = st.sidebar.text_area("Paste Raw DNA Sequence (A, C, G, T):", height=150)
    if pasted:
        engine = BruteForceRepetitionEngine(pasted, header="Pasted DNA Sequence")

# Sidebar: Brute Force Engine Parameters
st.sidebar.header("🔬 Brute-Force Parameters")

motif_range = st.sidebar.slider(
    "Combination Unit Length (Motif Size, k):",
    min_value=3,
    max_value=6,
    value=(3, 6),
    help="Generates all 4^k combinations for the selected lengths (3, 4, 5, 6 letters).",
)

rep_range = st.sidebar.slider(
    "Repetition Length Window (Total Letters):",
    min_value=6,
    max_value=10,
    value=(6, 10),
    help="Filters repetitions having total length within this window (6 to 10 letters).",
)

analysis_mode = st.sidebar.selectbox(
    "Repetition Mode:",
    ["Tandem Repeats (Microsatellites / SSRs)", "Interspersed Repeated Sequences", "Both"],
    index=2,
)

if engine:
    # Top Metrics Row
    m1, m2, m3, m4 = st.columns(4)
    with m1:
        st.metric("Sequence Length", f"{engine.length:,} nt")
    with m2:
        st.metric("GC Content", f"{engine.gc_content:.2f}%")
    with m3:
        # Calculate combinations generated
        combos_dict = engine.generate_combinations(motif_range[0], motif_range[1])
        total_combos = sum(len(c) for c in combos_dict.values())
        st.metric("Brute-Force Combinations", f"{total_combos:,}")
    with m4:
        st.metric("Window Size", f"{rep_range[0]} - {rep_range[1]} nt")

    st.caption(f"**Sequence Header:** `{engine.header}`")

    # Run Analysis
    with st.spinner("Executing Brute-Force Search Engine..."):
        tandem_repeats = []
        interspersed_repeats = []

        if analysis_mode in ["Tandem Repeats (Microsatellites / SSRs)", "Both"]:
            tandem_repeats = engine.detect_tandem_repeats(
                min_motif=motif_range[0],
                max_motif=motif_range[1],
                min_total_len=rep_range[0],
                max_total_len=rep_range[1],
            )

        if analysis_mode in ["Interspersed Repeated Sequences", "Both"]:
            interspersed_repeats = engine.detect_interspersed_repeats(
                min_len=rep_range[0],
                max_len=rep_range[1],
                min_occurrences=2,
            )

    # Tabs for Visualizations and Tables
    tab_overview, tab_tandem, tab_interspersed, tab_downloads = st.tabs([
        "📊 Visual Analytics",
        f"🔁 Tandem Repeats ({len(tandem_repeats)})",
        f"🔍 Interspersed Repeats ({len(interspersed_repeats)})",
        "💾 Data Export",
    ])

    with tab_overview:
        st.subheader("Spatial Distribution along Sequence")

        # Create Map Visualization
        fig, (ax_map, ax_hist) = plt.subplots(
            2, 1, figsize=(12, 7), facecolor="#0b0f19", gridspec_kw={"height_ratios": [1.4, 1.0]}
        )
        ax_map.set_facecolor("#0f172a")
        ax_hist.set_facecolor("#0f172a")

        color_map = {3: "#38bdf8", 4: "#34d399", 5: "#fbbf24", 6: "#f43f5e"}

        # Draw backbone
        ax_map.plot([1, engine.length], [0, 0], color="#475569", lw=4)

        for t in tandem_repeats:
            x = t["start_position"]
            y = t["total_length"]
            c = color_map.get(t["motif_length"], "#a855f7")
            ax_map.vlines(x, 0, y, color=c, alpha=0.6, lw=1.2)
            ax_map.scatter(x, y, color=c, s=30, edgecolors="#ffffff", lw=0.4)

        ax_map.set_xlim(-20, engine.length + 50)
        ax_map.set_ylim(-0.5, 12)
        ax_map.set_xlabel("Nucleotide Position (1 to N)", color="#cbd5e1", fontsize=9)
        ax_map.set_ylabel("Repetition Length (nt)", color="#cbd5e1", fontsize=9)
        ax_map.set_title("A. Genomic Positions of Detected Tandem Repetitions", color="#f8fafc", fontsize=11, fontweight="bold")
        ax_map.grid(True, ls=":", color="#334155", alpha=0.5)

        # Histogram
        len_counts = {}
        for l in range(rep_range[0], rep_range[1] + 1):
            len_counts[l] = sum(1 for t in tandem_repeats if t["total_length"] == l)

        all_l = list(range(rep_range[0], rep_range[1] + 1))
        counts = [len_counts[l] for l in all_l]

        bars = ax_hist.bar(all_l, counts, color="#6366f1", width=0.5, edgecolor="#a5b4fc")
        for bar, c in zip(bars, counts):
            if c > 0:
                ax_hist.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.8, str(c),
                             ha="center", va="bottom", color="#e0e7ff", fontsize=9, fontweight="bold")

        ax_hist.set_xticks(all_l)
        ax_hist.set_xticklabels([f"{l} letters" for l in all_l], color="#cbd5e1")
        ax_hist.set_xlabel("Repetition Length Window", color="#cbd5e1", fontsize=9)
        ax_hist.set_ylabel("Detected Count", color="#cbd5e1", fontsize=9)
        ax_hist.set_ylim(0, max(counts) + 10 if counts else 10)
        ax_hist.set_title("B. Repetition Frequency by Length", color="#f8fafc", fontsize=11, fontweight="bold")
        ax_hist.grid(axis="y", ls=":", color="#334155", alpha=0.5)

        plt.tight_layout()
        st.pyplot(fig)

    with tab_tandem:
        st.subheader("Tandem Repetition Results")
        if tandem_repeats:
            df_tandem = pd.DataFrame(tandem_repeats)
            st.dataframe(
                df_tandem[["pattern", "motif", "motif_length", "repeat_count", "total_length", "start_position", "end_position", "gc_content_pct"]],
                use_container_width=True,
                height=400,
            )
        else:
            st.info("No tandem repeats detected within the chosen parameters.")

    with tab_interspersed:
        st.subheader("Interspersed Repeated Sequences (Length 6–10 nt, Freq ≥ 2)")
        if interspersed_repeats:
            df_inter = pd.DataFrame([
                {
                    "pattern": r["pattern"],
                    "length": r["length"],
                    "occurrences_count": r["occurrences_count"],
                    "gc_content_pct": r["gc_content_pct"],
                    "first_position": r["positions"][0],
                    "all_positions": ", ".join(str(p) for p in r["positions"][:8]) + ("..." if len(r["positions"]) > 8 else ""),
                }
                for r in interspersed_repeats
            ])
            st.dataframe(df_inter, use_container_width=True, height=400)
        else:
            st.info("No interspersed repeats detected.")

    with tab_downloads:
        st.subheader("Export Results")
        c1, c2 = st.columns(2)
        with c1:
            if tandem_repeats:
                csv_data = pd.DataFrame(tandem_repeats).to_csv(index=False).encode("utf-8")
                st.download_button("📥 Download Tandem Repeats (CSV)", csv_data, "tandem_repeats.csv", "text/csv")
        with c2:
            if interspersed_repeats:
                csv_inter = pd.DataFrame(interspersed_repeats).to_csv(index=False).encode("utf-8")
                st.download_button("📥 Download Interspersed Repeats (CSV)", csv_inter, "interspersed_repeats.csv", "text/csv")
else:
    st.info("Please provide or load a DNA sequence from the sidebar to begin analysis.")
