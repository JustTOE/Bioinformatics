"""
app.py - Interactive Streamlit Web Application for Codon Analysis
Run with:
    streamlit run app.py
"""

import os
import sys
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

# Ensure local modules are found
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from genetic_code import codon_to_amino_acid, get_codon_info, translate_sequence, reverse_complement
from ncbi_search import CURATED_GENOMES, fetch_ncbi_fasta, search_ncbi_genomes
from codon_analyzer import analyze_frequencies, plot_top_10_codons, CLASS_COLORS

st.set_page_config(
    page_title="Codon Frequency Analyzer & Translator",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom Styling
st.markdown("""
<style>
    .main-title { font-size: 2.3rem; font-weight: 700; color: #1e3a8a; margin-bottom: 0.2rem; }
    .sub-title { font-size: 1.1rem; color: #475569; margin-bottom: 1.5rem; }
    .metric-box { background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 8px; padding: 12px; }
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-title">🧬 Codon Frequency & Genetic Code Analyzer</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-title">Multi-frame codon frequency detector, NCBI genome search, and genetic translation toolkit</div>', unsafe_allow_html=True)

# ----------------- SIDEBAR -----------------
st.sidebar.header("⚙️ Configuration & Data Source")

dataset_choice = st.sidebar.radio(
    "Choose Input Dataset:",
    [
        "Local E. coli Genomic DNA (Whole Chromosome)",
        "Local E. coli Coding Sequences (CDS)",
        "NCBI Curated Viral / Bacterial Genomes",
        "NCBI Accession / Live Search",
        "Upload Custom FASTA File",
        "Paste Raw DNA/RNA Sequence"
    ]
)

# Default path locations for the user's downloaded E. coli dataset
DEFAULT_ECOLI_GENOME = "/Users/toe/Documents/University/2026-2027/bioinfo/genome_datasets/ESCHERICHIA COLI/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
DEFAULT_ECOLI_CDS = "/Users/toe/Documents/University/2026-2027/bioinfo/genome_datasets/ESCHERICHIA COLI/ncbi_dataset/data/GCF_000005845.2/cds_from_genomic.fna"

temp_fasta_path = None
dataset_title = ""

if dataset_choice == "Local E. coli Genomic DNA (Whole Chromosome)":
    dataset_title = "Escherichia coli str. K-12 substr. MG1655 (Genomic DNA)"
    temp_fasta_path = DEFAULT_ECOLI_GENOME
    st.sidebar.info(f"Using local file:\n`{os.path.basename(DEFAULT_ECOLI_GENOME)}`")

elif dataset_choice == "Local E. coli Coding Sequences (CDS)":
    dataset_title = "Escherichia coli str. K-12 substr. MG1655 (CDS - 4,318 Genes)"
    temp_fasta_path = DEFAULT_ECOLI_CDS
    st.sidebar.info(f"Using local file:\n`{os.path.basename(DEFAULT_ECOLI_CDS)}`")

elif dataset_choice == "NCBI Curated Viral / Bacterial Genomes":
    genome_key = st.sidebar.selectbox(
        "Select Reference Genome:",
        list(CURATED_GENOMES.keys()),
        format_func=lambda k: f"{CURATED_GENOMES[k]['name']} ({CURATED_GENOMES[k]['accession']})"
    )
    g_info = CURATED_GENOMES[genome_key]
    dataset_title = f"{g_info['name']} ({g_info['accession']})"
    st.sidebar.caption(f"**Type:** {g_info['type']} | **Size:** {g_info['approx_size']}")
    st.sidebar.caption(g_info['description'])
    
    cache_dir = os.path.join(current_dir, "cache")
    os.makedirs(cache_dir, exist_ok=True)
    temp_fasta_path = os.path.join(cache_dir, f"{g_info['accession']}.fasta")
    if not os.path.exists(temp_fasta_path):
        with st.spinner(f"Fetching {g_info['accession']} from NCBI..."):
            fetch_ncbi_fasta(g_info['accession'], temp_fasta_path)

elif dataset_choice == "NCBI Accession / Live Search":
    acc_input = st.sidebar.text_input("NCBI Accession (e.g., NC_045512.2, NC_001416.1):", "NC_045512.2")
    if acc_input:
        dataset_title = f"NCBI Accession: {acc_input.strip()}"
        cache_dir = os.path.join(current_dir, "cache")
        os.makedirs(cache_dir, exist_ok=True)
        temp_fasta_path = os.path.join(cache_dir, f"{acc_input.strip()}.fasta")
        if not os.path.exists(temp_fasta_path):
            if st.sidebar.button("Fetch Genome from NCBI"):
                with st.spinner(f"Downloading {acc_input}..."):
                    try:
                        fetch_ncbi_fasta(acc_input.strip(), temp_fasta_path)
                        st.sidebar.success("Downloaded successfully!")
                    except Exception as e:
                        st.sidebar.error(f"Download failed: {e}")

elif dataset_choice == "Upload Custom FASTA File":
    uploaded = st.sidebar.file_uploader("Upload FASTA file (.fasta, .fna, .fa, .txt)", type=["fasta", "fna", "fa", "txt"])
    if uploaded is not None:
        cache_dir = os.path.join(current_dir, "cache")
        os.makedirs(cache_dir, exist_ok=True)
        temp_fasta_path = os.path.join(cache_dir, uploaded.name)
        with open(temp_fasta_path, "wb") as f:
            f.write(uploaded.getbuffer())
        dataset_title = f"Uploaded File: {uploaded.name}"

elif dataset_choice == "Paste Raw DNA/RNA Sequence":
    seq_input = st.sidebar.text_area("Paste nucleotide sequence:", "ATGAGCGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGC")
    if seq_input.strip():
        cache_dir = os.path.join(current_dir, "cache")
        os.makedirs(cache_dir, exist_ok=True)
        temp_fasta_path = os.path.join(cache_dir, "manual_input.fasta")
        with open(temp_fasta_path, "w") as f:
            f.write(f">Manual_Input_Sequence\n{seq_input.strip()}\n")
        dataset_title = "Manual Nucleotide Sequence"

st.sidebar.subheader("Reading Frame Selection")
frame_choice = st.sidebar.selectbox(
    "Analyze Frames:",
    ["All 6 Frames (+1, +2, +3, -1, -2, -3)", "Forward Frames (+1, +2, +3)", "Reverse Frames (-1, -2, -3)", "+1", "+2", "+3", "-1", "-2", "-3", "CDS Frame (+1 for each coding gene)"]
)

frame_mapping = {
    "All 6 Frames (+1, +2, +3, -1, -2, -3)": "all",
    "Forward Frames (+1, +2, +3)": "forward",
    "Reverse Frames (-1, -2, -3)": "reverse",
    "+1": "+1",
    "+2": "+2",
    "+3": "+3",
    "-1": "-1",
    "-2": "-2",
    "-3": "-3",
    "CDS Frame (+1 for each coding gene)": "cds"
}
frame_param = frame_mapping[frame_choice]

genetic_table = st.sidebar.selectbox("Genetic Code Translation Table:", [1, 11], format_func=lambda x: "Table 1 - Standard Code" if x == 1 else "Table 11 - Bacterial & Archaeal")

# ----------------- TABS -----------------
tab1, tab2, tab3, tab4 = st.tabs([
    "📊 Codon Frequency & Top 10 Chart",
    "🔤 Genetic Code Converter & Translator",
    "🔍 NCBI Search Tool",
    "🤖 LLM Prompt & Food Nutrition"
])

with tab1:
    if temp_fasta_path and os.path.exists(temp_fasta_path):
        st.subheader(f"Results for: {dataset_title}")
        
        with st.spinner("Calculating codon frequencies across requested frames..."):
            res = analyze_frequencies(temp_fasta_path, frame_mode=frame_param, table=genetic_table)
            
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Total Nucleotides", f"{res['total_nucleotides']:,} bp")
        c2.metric("Total Codons Counted", f"{res['total_codons']:,}")
        c3.metric("GC Content", f"{res['gc_content_percent']}%")
        top_codon = res['ranked_codons'][0]
        c4.metric("Top Codon", f"{top_codon['codon']} ({top_codon['amino_acid_3']})", f"{top_codon['frequency_percent']}%")
        
        # Display Chart
        st.markdown("### 🏆 Top 10 Most Frequent Codons")
        chart_file = os.path.join(current_dir, "cache", f"chart_{frame_param}.png")
        plot_top_10_codons(
            res['ranked_codons'],
            title=f"Top 10 Most Frequent Codons ({dataset_title.split('(')[0].strip()})",
            subtitle=f"Frame mode: {frame_choice} | Total Codons: {res['total_codons']:,} | GC: {res['gc_content_percent']}%",
            output_path=chart_file
        )
        st.image(chart_file, use_container_width=True)
        
        # Frequency Table
        st.markdown("### 📋 Complete Codon Frequency Table (Highest to Lowest)")
        df = pd.DataFrame(res['ranked_codons'])
        df = df.rename(columns={
            "rank": "Rank",
            "codon": "Codon (DNA)",
            "codon_rna": "Codon (RNA)",
            "amino_acid": "AA (1-Letter)",
            "amino_acid_3": "AA (3-Letter)",
            "amino_acid_name": "Amino Acid Name",
            "category": "Chemical Class",
            "count": "Total Count",
            "frequency_percent": "Frequency (%)",
            "frequency_per_thousand": "Per Thousand"
        })
        st.dataframe(df, use_container_width=True, height=400)
        
        csv_data = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download Full Codon Frequencies as CSV",
            data=csv_data,
            file_name=f"codon_frequencies_{frame_param}.csv",
            mime="text/csv"
        )
    else:
        st.warning("Please select or provide a FASTA file in the sidebar to begin analysis.")

with tab2:
    st.subheader("Interactive Codon to Amino Acid Converter")
    st.markdown("Enter any 3-nucleotide triplet to look up its corresponding amino acid, full chemical properties, and polarity.")
    
    col_input, col_out = st.columns([1, 2])
    with col_input:
        user_codon = st.text_input("Enter 3-letter Codon (DNA or RNA):", "ATG", max_chars=3)
        table_selected = st.radio("Translation Table:", [1, 11], format_func=lambda x: f"Table {x}")
        
    with col_out:
        if len(user_codon.strip()) == 3:
            info = get_codon_info(user_codon, table=table_selected)
            st.success(f"**Codon:** `{info['codon_dna']}` / `{info['codon_rna']}`")
            m1, m2, m3 = st.columns(3)
            m1.metric("Amino Acid", f"{info['amino_acid_name']} ({info['amino_acid_3']})")
            m2.metric("1-Letter Code", info['amino_acid_1'])
            m3.metric("Type", "Start Codon" if info['is_start'] else ("Stop Codon" if info['is_stop'] else "Sense Codon"))
            st.info(f"**Chemical Category:** {info['category']}")
        else:
            st.warning("Please enter a valid 3-character codon.")
            
    st.markdown("---")
    st.subheader("DNA/RNA Sequence Translation")
    test_seq = st.text_area("Enter Nucleotide Sequence to Translate:", "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA")
    trans_frame = st.selectbox("Reading Frame for Translation:", [+1, +2, +3, -1, -2, -3])
    if test_seq.strip():
        translated_prot = translate_sequence(test_seq, frame=trans_frame, table=table_selected)
        st.markdown(f"**Protein Sequence (Length: {len(translated_prot)} aa):**")
        st.code(translated_prot, language="text")

with tab3:
    st.subheader("NCBI Genome Search Tool")
    st.markdown("Search NCBI Nucleotide/RefSeq database directly for viral or bacterial genomes.")
    
    search_query = st.text_input("Enter search query:", "Bacteriophage T4 complete genome")
    if st.button("Search NCBI"):
        with st.spinner("Querying NCBI Entrez E-utilities..."):
            try:
                results = search_ncbi_genomes(search_query, max_results=6)
                if results:
                    st.success(f"Found {len(results)} records:")
                    for r in results:
                        st.markdown(f"**{r['accession']}**: {r['title']} *({int(r['length_bp']):,} bp)*")
                else:
                    st.info("No records found matching query.")
            except Exception as e:
                st.error(f"Search failed: {e}")

with tab4:
    st.subheader("🤖 LLM Prompt Generator & Amino Acid Food Finder")
    st.markdown(
        "Generate a structured prompt containing the **Top 10 Codons** and their translated amino acids, "
        "and determine the dietary food sources containing the **lowest frequency amino acid** among the Top 10."
    )
    
    prompt_mode = st.radio(
        "Choose Top 10 Context:",
        [
            "E. coli Whole Genome (All 6 Frames) -> Lowest Codon: GGC (Glycine)",
            "E. coli Coding Sequences (CDS) -> Lowest Codon: GTG (Valine)"
        ]
    )
    
    if "Glycine" in prompt_mode:
        prompt_file = os.path.join(current_dir, "prompt_whole_genome_glycine.txt")
        target_aa = "Glycine (Gly / G)"
        target_codon = "GGC (1.995%)"
        food_data = {
            "Animal-Based Foods": [
                ("Gelatin / Bone Broth Powder", "19,000 - 27,000 mg (19-27g)", "Highest known dietary source; dense in collagen peptides"),
                ("Pork Skin / Cracklings", "11,500 mg (11.5g)", "Rich in connective tissue collagen"),
                ("Chicken Skin / Roasted Wings", "3,300 mg (3.3g)", "Collagen-rich poultry cuts"),
                ("Beef Chuck / Brisket (Braised)", "2,100 mg (2.1g)", "Tougher cuts with high collagen breakdown"),
                ("Salmon / Whitefish with skin", "1,700 mg (1.7g)", "Fish collagen and lean marine protein"),
                ("Whole Eggs (Hard-boiled)", "680 mg (0.68g)", "Balanced essential and non-essential amino acids")
            ],
            "Plant-Based & Vegan Foods": [
                ("Spirulina (Dried Seaweed)", "3,100 mg (3.1g)", "Extremely dense microalgae protein"),
                ("Isolated Soy Protein", "3,600 mg (3.6g)", "Defatted soy protein isolate"),
                ("Sesame Seeds / Tahini", "2,150 mg (2.15g)", "Rich in minerals and sulfur-containing amino acids"),
                ("Pumpkin & Squash Seeds", "1,850 mg (1.85g)", "High zinc, magnesium, and glycine"),
                ("Hemp Seeds", "1,600 mg (1.6g)", "Complete plant protein profile"),
                ("Lentils & Chickpeas (Cooked)", "850 - 950 mg (0.85-0.95g)", "High-fiber legume staple")
            ]
        }
    else:
        prompt_file = os.path.join(current_dir, "prompt_cds_valine.txt")
        target_aa = "Valine (Val / V)"
        target_codon = "GTG (2.633%)"
        food_data = {
            "Animal-Based Foods": [
                ("Parmesan Cheese (Aged)", "2,050 mg (2.05g)", "High concentration of branched-chain amino acids (BCAAs)"),
                ("Beef Steak / Lean Sirloin", "1,850 mg (1.85g)", "Excellent BCAA density for muscle protein synthesis"),
                ("Chicken Breast (Roasted)", "1,700 mg (1.7g)", "High leucine, isoleucine, and valine triad"),
                ("Tuna / Wild Salmon", "1,550 mg (1.55g)", "Marine protein with omega-3 fatty acids"),
                ("Greek Yogurt / Cottage Cheese", "850 - 1,100 mg (0.85-1.1g)", "Slow-digesting casein rich in valine"),
                ("Whole Eggs", "860 mg (0.86g)", "Reference standard for amino acid bioavailability")
            ],
            "Plant-Based & Vegan Foods": [
                ("Soybeans / Edamame (Boiled)", "1,200 mg (1.2g)", "Complete plant protein with balanced BCAA profile"),
                ("Pumpkin & Watermelon Seeds", "1,550 mg (1.55g)", "Nutrient-dense seed source"),
                ("Nutritional Yeast", "1,800 mg (1.8g)", "Deactivated yeast fortified with B-vitamins"),
                ("Lentils & Black Beans (Cooked)", "650 - 750 mg (0.65-0.75g)", "Legume staple rich in lysine and valine"),
                ("Quinoa (Cooked)", "320 mg (0.32g)", "Complete grain containing all 9 essential amino acids"),
                ("Almonds & Walnuts", "700 - 850 mg (0.7-0.85g)", "Heart-healthy fats and micronutrients")
            ]
        }
        
    st.info(f"Target Amino Acid: **{target_aa}** | Associated with Codon: **{target_codon}**")
    
    if os.path.exists(prompt_file):
        with open(prompt_file, "r") as pf:
            prompt_text = pf.read()
            
        with st.expander("📄 View & Copy Full Generated LLM Prompt", expanded=False):
            st.code(prompt_text, language="markdown")
            
    st.markdown(f"### 🥗 Foods Rich in {target_aa}")
    col_a, col_p = st.columns(2)
    
    with col_a:
        st.markdown("#### 🥩 Animal-Based Sources")
        for food, amt, desc in food_data["Animal-Based Foods"]:
            st.markdown(f"- **{food}**: `{amt}` per 100g  \n  *{desc}*")
            
    with col_p:
        st.markdown("#### 🌱 Plant-Based & Vegan Sources")
        for food, amt, desc in food_data["Plant-Based & Vegan Foods"]:
            st.markdown(f"- **{food}**: `{amt}` per 100g  \n  *{desc}*")
            
    st.markdown("---")
    st.markdown("### 🔬 Biological Connection: Microbial Genetics & Nutrition")
    if "Glycine" in target_aa:
        st.markdown(
            "**In *E. coli*:** Glycine is the smallest amino acid with no chiral center. "
            "It confers structural flexibility in enzyme active sites, beta-turns, and the bacterial peptidoglycan cell wall (cross-linking glycan strands).  \n"
            "**In Humans:** Glycine is conditionally essential, crucial for collagen synthesis, glutathione (master antioxidant) production, bile acid conjugation, and central nervous system inhibitory neurotransmission."
        )
    else:
        st.markdown(
            "**In *E. coli*:** Valine is a branched-chain amino acid (BCAA) synthesized via acetolactate synthase. "
            "It is tightly regulated by feedback inhibition; its hydrophobic side chain is vital for packing bacterial core protein domains and transmembrane helices.  \n"
            "**In Humans:** Valine is one of the 9 essential amino acids (cannot be synthesized de novo). "
            "It stimulates muscle protein synthesis via mTOR signaling, assists in energy production during endurance exercise, and repairs damaged muscle tissue."
        )
