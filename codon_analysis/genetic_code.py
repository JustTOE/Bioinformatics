"""
genetic_code.py - Codon to Amino Acid Translation Module
Complies with NCBI Standard Genetic Code (Table 1) and Bacterial Code (Table 11).
"""

from typing import Dict, Tuple, Optional

# Standard Genetic Code (NCBI Translation Table 1)
STANDARD_GENETIC_CODE: Dict[str, str] = {
    # U / T
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    # C
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    # A
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    # G
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Bacterial, Archaeal and Plant Plastid Code (NCBI Translation Table 11)
# Note: In Table 11, the amino acid mapping is the same as Table 1,
# but alternative start codons include TTG, CTG, ATT, ATC, ATA, GTG.
BACTERIAL_GENETIC_CODE = STANDARD_GENETIC_CODE.copy()

# Amino Acid Metadata: (3-letter abbreviation, full name, chemical category)
AMINO_ACID_PROPERTIES = {
    'A': ('Ala', 'Alanine', 'Hydrophobic / Nonpolar'),
    'R': ('Arg', 'Arginine', 'Basic / Positively Charged'),
    'N': ('Asn', 'Asparagine', 'Polar / Uncharged'),
    'D': ('Asp', 'Aspartic Acid', 'Acidic / Negatively Charged'),
    'C': ('Cys', 'Cysteine', 'Polar / Special (Disulfide)'),
    'E': ('Glu', 'Glutamic Acid', 'Acidic / Negatively Charged'),
    'Q': ('Gln', 'Glutamine', 'Polar / Uncharged'),
    'G': ('Gly', 'Glycine', 'Special / Smallest'),
    'H': ('His', 'Histidine', 'Basic / Aromatic'),
    'I': ('Ile', 'Isoleucine', 'Hydrophobic / Nonpolar'),
    'L': ('Leu', 'Leucine', 'Hydrophobic / Nonpolar'),
    'K': ('Lys', 'Lysine', 'Basic / Positively Charged'),
    'M': ('Met', 'Methionine', 'Hydrophobic / Start'),
    'F': ('Phe', 'Phenylalanine', 'Hydrophobic / Aromatic'),
    'P': ('Pro', 'Proline', 'Special / Cyclic'),
    'S': ('Ser', 'Serine', 'Polar / Uncharged'),
    'T': ('Thr', 'Threonine', 'Polar / Uncharged'),
    'W': ('Trp', 'Tryptophan', 'Hydrophobic / Aromatic'),
    'Y': ('Tyr', 'Tyrosine', 'Polar / Aromatic'),
    'V': ('Val', 'Valine', 'Hydrophobic / Nonpolar'),
    '*': ('Stop', 'Translation Termination', 'Stop Codon'),
    'X': ('Unk', 'Unknown / Ambiguous', 'Ambiguous'),
}


def codon_to_amino_acid(codon: str, table: int = 1, strict: bool = False) -> str:
    """
    Translates a 3-nucleotide DNA or RNA codon into its 1-letter amino acid code.
    
    Parameters:
        codon (str): A 3-character string representing the codon (e.g., 'ATG', 'AUG', 'taa').
        table (int): NCBI translation table identifier (1 = Standard, 11 = Bacterial).
        strict (bool): If True, raises ValueError on invalid/ambiguous codons.
                       If False, returns 'X' for ambiguous codons.
                       
    Returns:
        str: 1-letter amino acid code (e.g. 'M') or '*' for stop codons.
        
    Examples:
        >>> codon_to_amino_acid('ATG')
        'M'
        >>> codon_to_amino_acid('aug')
        'M'
        >>> codon_to_amino_acid('TAA')
        '*'
    """
    clean_codon = codon.strip().upper().replace('U', 'T')
    
    if len(clean_codon) != 3:
        if strict:
            raise ValueError(f"Invalid codon length: '{codon}'. Expected exactly 3 nucleotides.")
        return 'X'
        
    code_table = BACTERIAL_GENETIC_CODE if table == 11 else STANDARD_GENETIC_CODE
    
    if clean_codon in code_table:
        return code_table[clean_codon]
    else:
        if strict:
            raise ValueError(f"Unknown codon: '{codon}' (contains ambiguous or non-standard bases).")
        return 'X'


def get_codon_info(codon: str, table: int = 1) -> Dict[str, str]:
    """
    Returns comprehensive biological information for a codon.
    """
    clean_codon = codon.strip().upper().replace('U', 'T')
    aa_1letter = codon_to_amino_acid(clean_codon, table=table, strict=False)
    abbr3, name, category = AMINO_ACID_PROPERTIES.get(aa_1letter, ('Unk', 'Unknown', 'Unknown'))
    
    return {
        'codon_dna': clean_codon,
        'codon_rna': clean_codon.replace('T', 'U'),
        'amino_acid_1': aa_1letter,
        'amino_acid_3': abbr3,
        'amino_acid_name': name,
        'category': category,
        'is_start': clean_codon in ('ATG', 'TTG', 'GTG') if table == 11 else clean_codon == 'ATG',
        'is_stop': aa_1letter == '*'
    }


def reverse_complement(sequence: str) -> str:
    """
    Calculates the 5' -> 3' reverse complement of a DNA sequence.
    """
    trans = str.maketrans('ATCGUNatcgun', 'TAGCNAntagcn')
    return sequence.translate(trans)[::-1]


def translate_sequence(sequence: str, frame: int = 1, table: int = 1) -> str:
    """
    Translates a DNA/RNA nucleotide sequence into a protein sequence in a given reading frame.
    
    Parameters:
        sequence (str): DNA/RNA sequence.
        frame (int): Reading frame: +1, +2, +3 for forward; -1, -2, -3 for reverse complement.
        table (int): NCBI genetic code table (1 = Standard, 11 = Bacterial).
        
    Returns:
        str: Protein sequence string.
    """
    seq = sequence.strip().upper().replace('U', 'T')
    
    if frame < 0:
        seq = reverse_complement(seq)
        offset = abs(frame) - 1
    else:
        offset = frame - 1
        
    if offset < 0 or offset >= 3:
        raise ValueError(f"Invalid frame: {frame}. Expected +1, +2, +3, -1, -2, or -3.")
        
    protein = []
    for i in range(offset, len(seq) - 2, 3):
        triplet = seq[i:i+3]
        protein.append(codon_to_amino_acid(triplet, table=table, strict=False))
        
    return ''.join(protein)
