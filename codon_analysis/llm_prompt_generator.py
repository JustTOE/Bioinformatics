"""
llm_prompt_generator.py - LLM Prompt Generator & Food Nutrition Query Engine
Generates an optimized LLM prompt from the Top 10 codon frequencies,
identifies the amino acid with the lowest frequency among the top 10,
and queries an LLM (online via Gemini/OpenAI REST API or locally via Ollama).
"""

import json
import os
import sys
import urllib.request
from typing import Dict, List, Optional, Tuple

# Import codon analysis utilities
try:
    from codon_analyzer import analyze_frequencies
    from genetic_code import get_codon_info
except ImportError:
    from .codon_analyzer import analyze_frequencies
    from .genetic_code import get_codon_info


def build_llm_prompt(ranked_codons: List[Dict], dataset_name: str = "Escherichia coli") -> Tuple[str, Dict]:
    """
    Constructs a structured, scientific prompt for an LLM containing the top 10 codons,
    their frequencies, and their translated amino acids, asking for the foods containing
    the lowest amino acid associated with the top 10 codons.
    """
    top_10 = ranked_codons[:10]
    
    # Identify the lowest frequency codon in the top 10
    lowest_codon_entry = top_10[-1]
    lowest_codon = lowest_codon_entry["codon"]
    lowest_aa_name = lowest_codon_entry["amino_acid_name"]
    lowest_aa_code = lowest_codon_entry["amino_acid"]
    lowest_aa_3 = lowest_codon_entry["amino_acid_3"]
    lowest_freq = lowest_codon_entry["frequency_percent"]
    
    # Format table for the prompt
    table_lines = [
        "| Rank | Codon (DNA) | Codon (RNA) | Amino Acid (1-Letter) | Amino Acid (3-Letter) | Full Name | Chemical Category | Frequency (%) | Total Counts |",
        "| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |"
    ]
    for row in top_10:
        table_lines.append(
            f"| {row['rank']} | {row['codon']} | {row['codon_rna']} | {row['amino_acid']} | "
            f"{row['amino_acid_3']} | {row['amino_acid_name']} | {row['category']} | "
            f"{row['frequency_percent']}% | {row['count']:,} |"
        )
    formatted_table = "\n".join(table_lines)
    
    prompt = f"""You are an expert computational biologist and clinical nutritionist.

### Context & Genomic Data:
We performed a genome-wide codon frequency analysis on **{dataset_name}**.
Below are the **Top 10 most frequent codons** detected, along with their respective translated amino acids, chemical classes, counts, and frequencies:

{formatted_table}

---

### Task & Objective:
1. **Identify the Lowest Amino Acid:**
   Based on the Top 10 codon table above:
   - Codon #{lowest_codon_entry['rank']} is `{lowest_codon}` (`{lowest_codon_entry['codon_rna']}`), which has the lowest frequency in this Top 10 group ({lowest_freq}%).
   - This codon codes for the amino acid: **{lowest_aa_name} ({lowest_aa_3} / {lowest_aa_code})**.

2. **Generate a Comprehensive Dietary Food List:**
   Provide a detailed, scientifically referenced list of foods that are rich sources of this specific amino acid (**{lowest_aa_name}**).
   
   Structure your answer with:
   - **Nutritional Profile of {lowest_aa_name}:** (Essential vs. non-essential / conditionally essential, recommended daily intake, biological role in human metabolism).
   - **Top Animal-Based Foods:** (Provide approximate content in mg or g per 100g serving).
   - **Top Plant-Based & Vegan Foods:** (Legumes, seeds, nuts, grains per 100g serving).
   - **Bioavailability & Absorption Factors:** (Digestibility, antagonistic amino acids, cooking impact).
   - **Biological Connection:** A brief concluding remark on why this amino acid is relevant both in bacterial protein expression ({dataset_name}) and human nutrition.
"""
    return prompt, lowest_codon_entry


def query_ollama_local(prompt: str, model: str = "llama3", base_url: str = "http://localhost:11434") -> Optional[str]:
    """
    Sends the prompt to a local Ollama instance (zero external dependency).
    """
    url = f"{base_url}/api/generate"
    payload = json.dumps({"model": model, "prompt": prompt, "stream": False}).encode("utf-8")
    req = urllib.request.Request(url, data=payload, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            return data.get("response")
    except Exception as e:
        return f"Local Ollama connection failed: {e}. (Ensure Ollama is running at {base_url})"


def query_gemini_api(prompt: str, api_key: Optional[str] = None) -> Optional[str]:
    """
    Sends the prompt to Google Gemini API via standard REST endpoint using urllib.
    """
    key = api_key or os.environ.get("GEMINI_API_KEY")
    if not key:
        return "GEMINI_API_KEY environment variable is not set."
        
    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key={key}"
    payload = json.dumps({
        "contents": [{"parts": [{"text": prompt}]}]
    }).encode("utf-8")
    
    req = urllib.request.Request(url, data=payload, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            return data["candidates"][0]["content"]["parts"][0]["text"]
    except Exception as e:
        return f"Gemini API error: {e}"


if __name__ == "__main__":
    # Test on E. coli whole genome
    ecoli_genome_path = "/Users/toe/Documents/University/2026-2027/bioinfo/genome_datasets/ESCHERICHIA COLI/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
    ecoli_cds_path = "/Users/toe/Documents/University/2026-2027/bioinfo/genome_datasets/ESCHERICHIA COLI/ncbi_dataset/data/GCF_000005845.2/cds_from_genomic.fna"
    
    print("Generating LLM Prompts for E. coli Top 10 Codons...")
    
    # 1. Whole Genome (6 frames)
    res_genome = analyze_frequencies(ecoli_genome_path, frame_mode="all")
    prompt_genome, lowest_genome = build_llm_prompt(res_genome["ranked_codons"], "Escherichia coli (Whole Genome 6-Frames)")
    
    out_file_genome = "/Users/toe/Documents/University/2026-2027/bioinfo/codon_analysis/prompt_whole_genome_glycine.txt"
    with open(out_file_genome, "w", encoding="utf-8") as f:
        f.write(prompt_genome)
    print(f"Saved Whole Genome prompt to: {out_file_genome}")
    print(f"Lowest Amino Acid in Top 10 (Whole Genome): {lowest_genome['amino_acid_name']} ({lowest_genome['codon']})")
    
    # 2. CDS (Coding Sequences)
    res_cds = analyze_frequencies(ecoli_cds_path, frame_mode="cds")
    prompt_cds, lowest_cds = build_llm_prompt(res_cds["ranked_codons"], "Escherichia coli (Coding Sequences / CDS)")
    
    out_file_cds = "/Users/toe/Documents/University/2026-2027/bioinfo/codon_analysis/prompt_cds_valine.txt"
    with open(out_file_cds, "w", encoding="utf-8") as f:
        f.write(prompt_cds)
    print(f"Saved CDS prompt to: {out_file_cds}")
    print(f"Lowest Amino Acid in Top 10 (CDS): {lowest_cds['amino_acid_name']} ({lowest_cds['codon']})")
