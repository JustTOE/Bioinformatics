"""
ncbi_search.py - NCBI Genome Search and Retrieval Module
Supports fetching viral and bacterial genomes via NCBI E-utilities (Entrez API)
and provides command-line guidance for NCBI Datasets CLI.
"""

import json
import os
import urllib.parse
import urllib.request
from typing import Dict, List, Optional, Tuple

# Curated reference genomes with RefSeq accessions
CURATED_GENOMES = {
    # Viral Genomes
    "sars-cov-2": {
        "name": "Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)",
        "type": "Viral",
        "accession": "NC_045512.2",
        "approx_size": "~29.9 kb",
        "description": "Positive-sense single-stranded RNA coronavirus"
    },
    "phage-lambda": {
        "name": "Enterobacteria phage lambda (Bacteriophage lambda)",
        "type": "Viral",
        "accession": "NC_001416.1",
        "approx_size": "~48.5 kb",
        "description": "Double-stranded DNA temperate phage infecting E. coli"
    },
    "phage-t4": {
        "name": "Escherichia phage T4 (Bacteriophage T4)",
        "type": "Viral",
        "accession": "NC_000866.4",
        "approx_size": "~168.9 kb",
        "description": "Large lytic dsDNA myovirus infecting E. coli"
    },
    "phix174": {
        "name": "Enterobacteria phage phiX174",
        "type": "Viral",
        "accession": "NC_001422.1",
        "approx_size": "~5.4 kb",
        "description": "Circular single-stranded DNA virus; first DNA genome ever sequenced"
    },
    # Bacterial Genomes
    "ecoli-mg1655": {
        "name": "Escherichia coli str. K-12 substr. MG1655",
        "type": "Bacterial",
        "accession": "NC_000913.3",
        "approx_size": "~4.64 Mb",
        "description": "Gold standard model bacterium for molecular biology and genomics"
    },
    "bacillus-subtilis": {
        "name": "Bacillus subtilis subsp. subtilis str. 168",
        "type": "Bacterial",
        "accession": "NC_000964.3",
        "approx_size": "~4.22 Mb",
        "description": "Model Gram-positive spore-forming bacterium"
    },
    "mycoplasma-genitalium": {
        "name": "Mycoplasma genitalium G37",
        "type": "Bacterial",
        "accession": "NC_000908.2",
        "approx_size": "~580 kb",
        "description": "One of the smallest known bacterial genomes"
    }
}


def search_ncbi_genomes(query: str, db: str = "nuccore", max_results: int = 5) -> List[Dict[str, str]]:
    """
    Searches NCBI for complete genomes matching a query.
    
    Parameters:
        query (str): Search term (e.g., 'SARS-CoV-2 complete genome' or 'Bacteriophage T4 complete genome').
        db (str): NCBI database to search ('nuccore' or 'assembly').
        max_results (int): Maximum number of records to return.
        
    Returns:
        List[Dict[str, str]]: List of search results with ID, title, and accession.
    """
    # Build ESearch URL
    search_term = f"{query} AND refseq[filter]"
    esearch_url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        f"db={db}&term={urllib.parse.quote(search_term)}&retmode=json&retmax={max_results}"
    )
    
    req = urllib.request.Request(esearch_url, headers={"User-Agent": "BioinfoCodonApp/1.0"})
    with urllib.request.urlopen(req, timeout=15) as resp:
        data = json.loads(resp.read().decode("utf-8"))
        
    id_list = data.get("esearchresult", {}).get("idlist", [])
    if not id_list:
        return []
        
    # Build ESummary URL to fetch titles and accessions
    ids_param = ",".join(id_list)
    esummary_url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        f"db={db}&id={ids_param}&retmode=json"
    )
    
    req_sum = urllib.request.Request(esummary_url, headers={"User-Agent": "BioinfoCodonApp/1.0"})
    with urllib.request.urlopen(req_sum, timeout=15) as resp_sum:
        summary_data = json.loads(resp_sum.read().decode("utf-8"))
        
    result_dict = summary_data.get("result", {})
    results = []
    for uid in id_list:
        item = result_dict.get(uid, {})
        caption = item.get("caption", uid)
        title = item.get("title", "No title available")
        slen = item.get("slen", "Unknown")
        results.append({
            "uid": uid,
            "accession": caption,
            "title": title,
            "length_bp": str(slen)
        })
        
    return results


def fetch_ncbi_fasta(accession: str, output_file: Optional[str] = None) -> Tuple[str, str]:
    """
    Downloads a complete nucleotide FASTA sequence from NCBI given a RefSeq/GenBank accession.
    
    Parameters:
        accession (str): NCBI accession number (e.g. 'NC_045512.2' or 'NC_000913.3').
        output_file (str, optional): If provided, writes FASTA to this path.
        
    Returns:
        Tuple[str, str]: (FASTA header, clean nucleotide sequence)
    """
    clean_acc = accession.strip()
    efetch_url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        f"db=nuccore&id={clean_acc}&rettype=fasta&retmode=text"
    )
    
    req = urllib.request.Request(efetch_url, headers={"User-Agent": "BioinfoCodonApp/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        fasta_text = resp.read().decode("utf-8")
        
    lines = fasta_text.strip().splitlines()
    header = lines[0] if lines and lines[0].startswith(">") else f">{clean_acc}"
    sequence = "".join([l.strip().upper() for l in lines if not l.startswith(">")])
    
    if output_file:
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        with open(output_file, "w") as f:
            f.write(f"{header}\n")
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")
                
    return header, sequence


if __name__ == "__main__":
    print("Testing NCBI search for SARS-CoV-2...")
    results = search_ncbi_genomes("SARS-CoV-2 complete genome")
    for r in results:
        print(f"[{r['accession']}] {r['title']} ({r['length_bp']} bp)")
