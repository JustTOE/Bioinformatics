"""
simulator.py - Artificial DNA Sequence Generator with Transposable Element Simulation
======================================================================================
Generates synthetic DNA sequences of 200 to 500 base pairs containing 3 to 4
transposable elements (TEs), with support for nested/intersected configurations.

Biologically modeled features:
- Class II DNA transposons with Terminal Inverted Repeats (TIRs) and Target Site Duplications (TSDs).
- Class I Retrotransposons with Long Terminal Repeats (LTRs) and TSDs.
- Nested transposable element insertions (transposon jumping into another transposon).
- Overlapping / chimeric transposon insertions.
- Complete ground truth tracking (exact coordinates, fragments, hierarchy).
"""

import random
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    return ''.join(complement.get(b, 'N') for b in reversed(seq))


@dataclass
class TransposonModel:
    """Definition of a canonical transposable element family."""
    name: str
    family_type: str            # 'DNA_TIR' or 'RETRO_LTR'
    core_length: int            # Length of internal sequence
    tir_or_ltr_seq: str         # 5' terminal sequence (TIR or LTR)
    tsd_motif: str              # Canonical target site duplication motif (e.g., 'TA', 'TTAA')
    color: str                  # Hex color for visualization
    full_sequence: str = ""     # Assembled full canonical reference sequence
    description: str = ""

    def __post_init__(self):
        if not self.full_sequence:
            self.full_sequence = self.generate_canonical_sequence()

    def generate_canonical_sequence(self, seed: Optional[int] = None) -> str:
        """Construct the canonical full-length sequence of this TE family."""
        rng = random.Random(seed if seed is not None else hash(self.name))
        
        # Internal core body
        core = "".join(rng.choices(["A", "C", "G", "T"], k=self.core_length))
        
        if self.family_type == "DNA_TIR":
            # 5' TIR + Core + 3' TIR (reverse complement of 5' TIR)
            tir_5 = self.tir_or_ltr_seq
            tir_3 = reverse_complement(tir_5)
            return tir_5 + core + tir_3
        elif self.family_type == "RETRO_LTR":
            # 5' LTR + Core + 3' LTR (direct repeat of 5' LTR)
            ltr_5 = self.tir_or_ltr_seq
            ltr_3 = ltr_5
            return ltr_5 + core + ltr_3
        else:
            return core


@dataclass
class GroundTruthFragment:
    """Represents a physically continuous fragment of a transposable element."""
    element_name: str
    start: int                  # 0-based inclusive start index in final genomic sequence
    end: int                    # 0-based inclusive end index in final genomic sequence
    ref_start: int              # 0-based start index within canonical reference sequence
    ref_end: int                # 0-based end index within canonical reference sequence
    fragment_label: str         # e.g., 'left_flank', 'right_flank', 'intact'
    sequence: str


@dataclass
class GroundTruthElement:
    """Represents a simulated transposable element insertion event."""
    element_id: str             # e.g., 'TE1_Mariner'
    family_name: str
    full_start: int             # 0-based overall outer start index
    full_end: int               # 0-based overall outer end index (inclusive)
    is_interrupted: bool = False
    interrupted_by: List[str] = field(default_factory=list)
    parent_element: Optional[str] = None
    nesting_level: int = 0      # 0 = root, 1 = nested once, 2 = nested twice
    tsd_sequence: str = ""
    tsd_start_5: Optional[int] = None
    tsd_start_3: Optional[int] = None
    fragments: List[GroundTruthFragment] = field(default_factory=list)
    canonical_sequence: str = ""


@dataclass
class SimulationResult:
    """Holds the complete output of an artificial sequence simulation."""
    dna_sequence: str
    total_length: int
    scenario_name: str
    scenario_description: str
    elements: Dict[str, GroundTruthElement]
    reference_library: Dict[str, str]       # name -> canonical sequence
    insertion_history: List[str]            # Log of insertion events


def get_default_te_library() -> Dict[str, TransposonModel]:
    """Define standard reference library of 4 realistic synthetic transposable elements."""
    return {
        "TE1_Mariner": TransposonModel(
            name="TE1_Mariner",
            family_type="DNA_TIR",
            core_length=50,
            tir_or_ltr_seq="CAGTGTCTC",     # 9 bp TIR; 3' is GAGACACTG
            tsd_motif="TA",                 # 2 bp TSD (classic Tc1/mariner)
            color="#2563EB",                # Blue
            description="Tc1/Mariner-like DNA transposon with 9-bp TIRs and 2-bp TA TSD"
        ),
        "TE2_PiggyBac": TransposonModel(
            name="TE2_PiggyBac",
            family_type="DNA_TIR",
            core_length=40,
            tir_or_ltr_seq="CCCTAGTT",      # 8 bp TIR; 3' is AACTAGGG
            tsd_motif="TTAA",               # 4 bp TSD (classic piggyBac)
            color="#DC2626",                # Red
            description="PiggyBac-like DNA transposon with 8-bp TIRs and 4-bp TTAA TSD"
        ),
        "TE3_Gypsy": TransposonModel(
            name="TE3_Gypsy",
            family_type="RETRO_LTR",
            core_length=42,
            tir_or_ltr_seq="TGTTGACC",      # 8 bp LTR direct repeat
            tsd_motif="AAGT",               # 4 bp TSD
            color="#16A34A",                # Green
            description="Gypsy/Ty3-like LTR retrotransposon with direct LTRs and 4-bp TSD"
        ),
        "TE4_Helitron": TransposonModel(
            name="TE4_Helitron",
            family_type="DNA_TIR",
            core_length=36,
            tir_or_ltr_seq="GCTCGTAG",      # 8 bp TIR
            tsd_motif="CGT",                # 3 bp TSD
            color="#9333EA",                # Purple
            description="Helitron/hAT-like transposable element with 8-bp TIRs and 3-bp TSD"
        )
    }


def generate_random_dna(length: int, gc_content: float = 0.50, seed: Optional[int] = None) -> str:
    """Generate a random DNA background with specified GC content."""
    rng = random.Random(seed)
    p_gc = gc_content / 2.0
    p_at = (1.0 - gc_content) / 2.0
    weights = [p_at, p_gc, p_gc, p_at]  # A, C, G, T
    return "".join(rng.choices(["A", "C", "G", "T"], weights=weights, k=length))


class TransposonSimulator:
    """
    Simulates artificial DNA sequences (200-500 bp) with 3 to 4 transposable elements,
    supporting nested (intersected) insertions and overlapping coordinates.
    """

    def __init__(self, te_library: Optional[Dict[str, TransposonModel]] = None, seed: Optional[int] = None):
        self.te_library = te_library or get_default_te_library()
        self.seed = seed
        self.rng = random.Random(seed)

    def simulate(self,
                 target_length: int = 350,
                 num_elements: int = 4,
                 intersection_type: str = "nested_2",
                 gc_content: float = 0.50) -> SimulationResult:
        """
        Simulate an artificial DNA sequence with 3 or 4 transposable elements.

        Parameters:
        -----------
        target_length: int, 200 to 500 bp
        num_elements: int, 3 or 4
        intersection_type:
            - 'nested_2': 2 transposons intersect (TE2 nested inside TE1)
            - 'nested_3': 3 transposons intersect (TE3 inside TE2 inside TE1 - Russian doll)
            - 'overlapping': 2 transposons partially overlap in coordinates
            - 'independent': 3 or 4 non-intersecting transposons (baseline control)
        gc_content: float, background GC proportion
        """
        assert 200 <= target_length <= 500, f"Target length {target_length} must be between 200 and 500 bp"
        assert num_elements in (3, 4), f"num_elements {num_elements} must be 3 or 4"

        # Select elements
        element_keys = ["TE1_Mariner", "TE2_PiggyBac", "TE3_Gypsy"]
        if num_elements == 4:
            element_keys.append("TE4_Helitron")

        selected_models = {k: self.te_library[k] for k in element_keys}

        if intersection_type == "nested_2":
            return self._simulate_nested_2(target_length, selected_models, gc_content)
        elif intersection_type == "nested_3":
            return self._simulate_nested_3(target_length, selected_models, gc_content)
        elif intersection_type == "overlapping":
            return self._simulate_overlapping(target_length, selected_models, gc_content)
        elif intersection_type == "independent":
            return self._simulate_independent(target_length, selected_models, gc_content)
        else:
            raise ValueError(f"Unknown intersection_type: {intersection_type}")

    def _simulate_nested_2(self,
                           target_length: int,
                           models: Dict[str, TransposonModel],
                           gc_content: float) -> SimulationResult:
        """
        Simulate 2 intersected transposons:
        TE2 inserts INSIDE TE1, interrupting it into left and right pieces.
        TE3 (and TE4 if requested) are placed downstream or upstream.
        """
        history = []
        te1 = models["TE1_Mariner"]
        te2 = models["TE2_PiggyBac"]
        te3 = models["TE3_Gypsy"]
        te4 = models.get("TE4_Helitron", None)

        seq_te1 = te1.full_sequence
        seq_te2 = te2.full_sequence
        seq_te3 = te3.full_sequence
        seq_te4 = te4.full_sequence if te4 else ""

        # Choose interruption point in TE1 (approximately in the middle, within the core)
        # TE1 structure: TIR_5 (9 bp) + core (50 bp) + TIR_3 (9 bp) = 68 bp total
        te1_split = len(te1.tir_or_ltr_seq) + self.rng.randint(18, 32)
        te1_left = seq_te1[:te1_split]
        te1_right = seq_te1[te1_split:]

        # Compound nested structure: te1_left + te2.tsd_motif + seq_te2 + te2.tsd_motif + te1_right
        # In biological insertion, target site duplication (TSD) is duplicated at the insertion boundary:
        tsd2 = te2.tsd_motif
        compound_te1_te2 = te1_left + tsd2 + seq_te2 + tsd2 + te1_right
        len_compound = len(compound_te1_te2)

        # Remaining elements length
        len_other_tes = len(seq_te3) + len(te3.tsd_motif) * 2
        if te4:
            len_other_tes += len(seq_te4) + len(te4.tsd_motif) * 2

        # Background host length to reach target_length
        total_te_len = len_compound + len_other_tes
        host_needed = max(target_length - total_te_len, 40)
        
        # Distribute host into 3 or 4 spacers:
        # [Host1] [Compound TE1/TE2] [Host2] [TE3] [Host3] [TE4 (optional)] [Host4]
        num_spacers = 4 if te4 else 3
        spacers = self._partition_integer(host_needed, num_spacers, min_val=8)

        # Assemble sequence with biological TSDs
        tsd1 = te1.tsd_motif
        tsd3 = te3.tsd_motif
        tsd4 = te4.tsd_motif if te4 else ""

        seq_parts = []
        # Host 1
        h1 = generate_random_dna(spacers[0], gc_content, self.rng.randint(1, 100000))
        seq_parts.append(h1)

        # Compound TE1/TE2 insertion
        # Outer TE1 has its TSD flanking it
        pos_tsd1_5 = len("".join(seq_parts))
        seq_parts.append(tsd1)
        start_te1 = len("".join(seq_parts))
        
        start_te1_left = start_te1
        seq_parts.append(te1_left)
        end_te1_left = len("".join(seq_parts)) - 1

        # TE2 insertion inside TE1
        pos_tsd2_5 = len("".join(seq_parts))
        seq_parts.append(tsd2)
        start_te2 = len("".join(seq_parts))
        seq_parts.append(seq_te2)
        end_te2 = len("".join(seq_parts)) - 1
        pos_tsd2_3 = len("".join(seq_parts))
        seq_parts.append(tsd2)

        start_te1_right = len("".join(seq_parts))
        seq_parts.append(te1_right)
        end_te1_right = len("".join(seq_parts)) - 1
        end_te1 = end_te1_right

        pos_tsd1_3 = len("".join(seq_parts))
        seq_parts.append(tsd1)

        # Host 2
        h2 = generate_random_dna(spacers[1], gc_content, self.rng.randint(1, 100000))
        seq_parts.append(h2)

        # TE3 insertion
        pos_tsd3_5 = len("".join(seq_parts))
        seq_parts.append(tsd3)
        start_te3 = len("".join(seq_parts))
        seq_parts.append(seq_te3)
        end_te3 = len("".join(seq_parts)) - 1
        pos_tsd3_3 = len("".join(seq_parts))
        seq_parts.append(tsd3)

        # Host 3
        h3 = generate_random_dna(spacers[2], gc_content, self.rng.randint(1, 100000))
        seq_parts.append(h3)

        if te4:
            pos_tsd4_5 = len("".join(seq_parts))
            seq_parts.append(tsd4)
            start_te4 = len("".join(seq_parts))
            seq_parts.append(seq_te4)
            end_te4 = len("".join(seq_parts)) - 1
            pos_tsd4_3 = len("".join(seq_parts))
            seq_parts.append(tsd4)

            # Host 4
            h4 = generate_random_dna(spacers[3], gc_content, self.rng.randint(1, 100000))
            seq_parts.append(h4)

        final_seq = "".join(seq_parts)

        # Build Ground Truth data structures
        elements = {}
        
        # TE1
        frag1_left = GroundTruthFragment(
            element_name="TE1_Mariner",
            start=start_te1_left,
            end=end_te1_left,
            ref_start=0,
            ref_end=te1_split - 1,
            fragment_label="5'_fragment",
            sequence=te1_left
        )
        frag1_right = GroundTruthFragment(
            element_name="TE1_Mariner",
            start=start_te1_right,
            end=end_te1_right,
            ref_start=te1_split,
            ref_end=len(seq_te1) - 1,
            fragment_label="3'_fragment",
            sequence=te1_right
        )
        elements["TE1_Mariner"] = GroundTruthElement(
            element_id="TE1_Mariner",
            family_name="TE1_Mariner",
            full_start=start_te1,
            full_end=end_te1,
            is_interrupted=True,
            interrupted_by=["TE2_PiggyBac"],
            parent_element=None,
            nesting_level=0,
            tsd_sequence=tsd1,
            tsd_start_5=pos_tsd1_5,
            tsd_start_3=pos_tsd1_3,
            fragments=[frag1_left, frag1_right],
            canonical_sequence=seq_te1
        )

        # TE2
        frag2 = GroundTruthFragment(
            element_name="TE2_PiggyBac",
            start=start_te2,
            end=end_te2,
            ref_start=0,
            ref_end=len(seq_te2) - 1,
            fragment_label="intact_nested",
            sequence=seq_te2
        )
        elements["TE2_PiggyBac"] = GroundTruthElement(
            element_id="TE2_PiggyBac",
            family_name="TE2_PiggyBac",
            full_start=start_te2,
            full_end=end_te2,
            is_interrupted=False,
            interrupted_by=[],
            parent_element="TE1_Mariner",
            nesting_level=1,
            tsd_sequence=tsd2,
            tsd_start_5=pos_tsd2_5,
            tsd_start_3=pos_tsd2_3,
            fragments=[frag2],
            canonical_sequence=seq_te2
        )

        # TE3
        frag3 = GroundTruthFragment(
            element_name="TE3_Gypsy",
            start=start_te3,
            end=end_te3,
            ref_start=0,
            ref_end=len(seq_te3) - 1,
            fragment_label="intact_independent",
            sequence=seq_te3
        )
        elements["TE3_Gypsy"] = GroundTruthElement(
            element_id="TE3_Gypsy",
            family_name="TE3_Gypsy",
            full_start=start_te3,
            full_end=end_te3,
            is_interrupted=False,
            interrupted_by=[],
            parent_element=None,
            nesting_level=0,
            tsd_sequence=tsd3,
            tsd_start_5=pos_tsd3_5,
            tsd_start_3=pos_tsd3_3,
            fragments=[frag3],
            canonical_sequence=seq_te3
        )

        # TE4 (if present)
        if te4:
            frag4 = GroundTruthFragment(
                element_name="TE4_Helitron",
                start=start_te4,
                end=end_te4,
                ref_start=0,
                ref_end=len(seq_te4) - 1,
                fragment_label="intact_independent",
                sequence=seq_te4
            )
            elements["TE4_Helitron"] = GroundTruthElement(
                element_id="TE4_Helitron",
                family_name="TE4_Helitron",
                full_start=start_te4,
                full_end=end_te4,
                is_interrupted=False,
                interrupted_by=[],
                parent_element=None,
                nesting_level=0,
                tsd_sequence=tsd4,
                tsd_start_5=pos_tsd4_5,
                tsd_start_3=pos_tsd4_3,
                fragments=[frag4],
                canonical_sequence=seq_te4
            )

        history.append("1. TE1_Mariner inserted into genomic background.")
        history.append(f"2. TE2_PiggyBac inserted INSIDE TE1_Mariner at relative base {te1_split}, splitting TE1 into 2 fragments.")
        history.append("3. TE3_Gypsy inserted independently downstream.")
        if te4:
            history.append("4. TE4_Helitron inserted independently downstream.")

        return SimulationResult(
            dna_sequence=final_seq,
            total_length=len(final_seq),
            scenario_name="2 Intersected Transposons (Single Nesting / Interruption)",
            scenario_description=(
                f"Simulates TE2_PiggyBac inserted directly inside TE1_Mariner (nested insertion). "
                f"TE1 is broken into two fragments ({len(te1_left)} bp 5' and {len(te1_right)} bp 3') "
                f"separated by TE2 ({len(seq_te2)} bp). Remaining element(s) are independent."
            ),
            elements=elements,
            reference_library={k: m.full_sequence for k, m in models.items()},
            insertion_history=history
        )

    def _simulate_nested_3(self,
                           target_length: int,
                           models: Dict[str, TransposonModel],
                           gc_content: float) -> SimulationResult:
        """
        Simulate 3 intersected transposons (Double Nesting / Russian Doll insertion):
        TE3 inserts INSIDE TE2, which itself is inserted INSIDE TE1!
        TE1 is split into 2 pieces; TE2 is split into 2 pieces; TE3 is intact at the core.
        TE4 (if present) is placed downstream.
        """
        history = []
        te1 = models["TE1_Mariner"]
        te2 = models["TE2_PiggyBac"]
        te3 = models["TE3_Gypsy"]
        te4 = models.get("TE4_Helitron", None)

        seq_te1 = te1.full_sequence
        seq_te2 = te2.full_sequence
        seq_te3 = te3.full_sequence
        seq_te4 = te4.full_sequence if te4 else ""

        # Split TE2 for TE3 insertion
        te2_split = len(te2.tir_or_ltr_seq) + self.rng.randint(12, 22)
        te2_left = seq_te2[:te2_split]
        te2_right = seq_te2[te2_split:]

        # Split TE1 for TE2 compound insertion
        te1_split = len(te1.tir_or_ltr_seq) + self.rng.randint(16, 26)
        te1_left = seq_te1[:te1_split]
        te1_right = seq_te1[te1_split:]

        tsd3 = te3.tsd_motif
        tsd2 = te2.tsd_motif
        tsd1 = te1.tsd_motif

        # Calculate space
        compound_inner_len = len(te2_left) + len(tsd3)*2 + len(seq_te3) + len(te2_right)
        compound_total_len = len(te1_left) + len(tsd2)*2 + compound_inner_len + len(te1_right)
        len_te4 = (len(seq_te4) + len(te4.tsd_motif)*2) if te4 else 0
        total_te_len = compound_total_len + len_te4

        host_needed = max(target_length - total_te_len, 35)
        num_spacers = 3 if te4 else 2
        spacers = self._partition_integer(host_needed, num_spacers, min_val=8)

        seq_parts = []

        # Host 1
        h1 = generate_random_dna(spacers[0], gc_content, self.rng.randint(1, 100000))
        seq_parts.append(h1)

        # Compound Nesting assembly:
        # TE1 Outer Left
        pos_tsd1_5 = len("".join(seq_parts))
        seq_parts.append(tsd1)
        start_te1 = len("".join(seq_parts))
        start_te1_left = start_te1
        seq_parts.append(te1_left)
        end_te1_left = len("".join(seq_parts)) - 1

        # TE2 Outer Left (nested inside TE1)
        pos_tsd2_5 = len("".join(seq_parts))
        seq_parts.append(tsd2)
        start_te2 = len("".join(seq_parts))
        start_te2_left = start_te2
        seq_parts.append(te2_left)
        end_te2_left = len("".join(seq_parts)) - 1

        # TE3 Innermost (nested inside TE2)
        pos_tsd3_5 = len("".join(seq_parts))
        seq_parts.append(tsd3)
        start_te3 = len("".join(seq_parts))
        seq_parts.append(seq_te3)
        end_te3 = len("".join(seq_parts)) - 1
        pos_tsd3_3 = len("".join(seq_parts))
        seq_parts.append(tsd3)

        # TE2 Outer Right
        start_te2_right = len("".join(seq_parts))
        seq_parts.append(te2_right)
        end_te2_right = len("".join(seq_parts)) - 1
        end_te2 = end_te2_right
        pos_tsd2_3 = len("".join(seq_parts))
        seq_parts.append(tsd2)

        # TE1 Outer Right
        start_te1_right = len("".join(seq_parts))
        seq_parts.append(te1_right)
        end_te1_right = len("".join(seq_parts)) - 1
        end_te1 = end_te1_right
        pos_tsd1_3 = len("".join(seq_parts))
        seq_parts.append(tsd1)

        # Host 2
        h2 = generate_random_dna(spacers[1], gc_content, self.rng.randint(1, 100000))
        seq_parts.append(h2)

        if te4:
            tsd4 = te4.tsd_motif
            pos_tsd4_5 = len("".join(seq_parts))
            seq_parts.append(tsd4)
            start_te4 = len("".join(seq_parts))
            seq_parts.append(seq_te4)
            end_te4 = len("".join(seq_parts)) - 1
            pos_tsd4_3 = len("".join(seq_parts))
            seq_parts.append(tsd4)

            # Host 3
            h3 = generate_random_dna(spacers[2], gc_content, self.rng.randint(1, 100000))
            seq_parts.append(h3)

        final_seq = "".join(seq_parts)

        # Ground truth
        elements = {}
        # TE1 (Nesting level 0, interrupted by TE2)
        elements["TE1_Mariner"] = GroundTruthElement(
            element_id="TE1_Mariner",
            family_name="TE1_Mariner",
            full_start=start_te1,
            full_end=end_te1,
            is_interrupted=True,
            interrupted_by=["TE2_PiggyBac"],
            parent_element=None,
            nesting_level=0,
            tsd_sequence=tsd1,
            tsd_start_5=pos_tsd1_5,
            tsd_start_3=pos_tsd1_3,
            fragments=[
                GroundTruthFragment("TE1_Mariner", start_te1_left, end_te1_left, 0, te1_split - 1, "5'_fragment", te1_left),
                GroundTruthFragment("TE1_Mariner", start_te1_right, end_te1_right, te1_split, len(seq_te1) - 1, "3'_fragment", te1_right)
            ],
            canonical_sequence=seq_te1
        )

        # TE2 (Nesting level 1, interrupted by TE3)
        elements["TE2_PiggyBac"] = GroundTruthElement(
            element_id="TE2_PiggyBac",
            family_name="TE2_PiggyBac",
            full_start=start_te2,
            full_end=end_te2,
            is_interrupted=True,
            interrupted_by=["TE3_Gypsy"],
            parent_element="TE1_Mariner",
            nesting_level=1,
            tsd_sequence=tsd2,
            tsd_start_5=pos_tsd2_5,
            tsd_start_3=pos_tsd2_3,
            fragments=[
                GroundTruthFragment("TE2_PiggyBac", start_te2_left, end_te2_left, 0, te2_split - 1, "5'_fragment", te2_left),
                GroundTruthFragment("TE2_PiggyBac", start_te2_right, end_te2_right, te2_split, len(seq_te2) - 1, "3'_fragment", te2_right)
            ],
            canonical_sequence=seq_te2
        )

        # TE3 (Nesting level 2, intact innermost core)
        elements["TE3_Gypsy"] = GroundTruthElement(
            element_id="TE3_Gypsy",
            family_name="TE3_Gypsy",
            full_start=start_te3,
            full_end=end_te3,
            is_interrupted=False,
            interrupted_by=[],
            parent_element="TE2_PiggyBac",
            nesting_level=2,
            tsd_sequence=tsd3,
            tsd_start_5=pos_tsd3_5,
            tsd_start_3=pos_tsd3_3,
            fragments=[
                GroundTruthFragment("TE3_Gypsy", start_te3, end_te3, 0, len(seq_te3) - 1, "intact_core", seq_te3)
            ],
            canonical_sequence=seq_te3
        )

        if te4:
            elements["TE4_Helitron"] = GroundTruthElement(
                element_id="TE4_Helitron",
                family_name="TE4_Helitron",
                full_start=start_te4,
                full_end=end_te4,
                is_interrupted=False,
                interrupted_by=[],
                parent_element=None,
                nesting_level=0,
                tsd_sequence=tsd4,
                tsd_start_5=pos_tsd4_5,
                tsd_start_3=pos_tsd4_3,
                fragments=[
                    GroundTruthFragment("TE4_Helitron", start_te4, end_te4, 0, len(seq_te4) - 1, "intact_independent", seq_te4)
                ],
                canonical_sequence=seq_te4
            )

        history.append("1. TE1_Mariner inserted into genomic background.")
        history.append(f"2. TE2_PiggyBac inserted INSIDE TE1_Mariner (Level 1 nesting).")
        history.append(f"3. TE3_Gypsy inserted INSIDE TE2_PiggyBac (Level 2 nested Russian doll).")
        if te4:
            history.append("4. TE4_Helitron inserted independently downstream.")

        return SimulationResult(
            dna_sequence=final_seq,
            total_length=len(final_seq),
            scenario_name="3 Intersected Transposons (Double Russian-Doll Nesting)",
            scenario_description=(
                f"Simulates 3 deeply intersected transposons: TE3_Gypsy is nested INSIDE TE2_PiggyBac, "
                f"which is nested INSIDE TE1_Mariner. TE1 and TE2 are both interrupted into two fragments, "
                f"while TE3 sits at the center."
            ),
            elements=elements,
            reference_library={k: m.full_sequence for k, m in models.items()},
            insertion_history=history
        )

    def _simulate_overlapping(self,
                             target_length: int,
                             models: Dict[str, TransposonModel],
                             gc_content: float) -> SimulationResult:
        """
        Simulate coordinate-overlapping transposons:
        TE1 and TE2 share a 15-25 bp overlapping region due to partial target deletion or
        tandem insertion recombination.
        """
        history = []
        te1 = models["TE1_Mariner"]
        te2 = models["TE2_PiggyBac"]
        te3 = models["TE3_Gypsy"]
        te4 = models.get("TE4_Helitron", None)

        seq_te1 = te1.full_sequence
        seq_te2 = te2.full_sequence
        seq_te3 = te3.full_sequence
        seq_te4 = te4.full_sequence if te4 else ""

        overlap_len = 18
        # Create an overlap where TE1 suffix and TE2 prefix share sequence or overlap in sequence
        # We merge TE1 and TE2 with an overlap of overlap_len bases
        overlap_seq = seq_te1[-overlap_len:]
        modified_te2 = overlap_seq + seq_te2[overlap_len:]

        combined_overlap = seq_te1 + seq_te2[overlap_len:]
        len_combined = len(combined_overlap)

        len_others = len(seq_te3) + (len(seq_te4) if te4 else 0)
        host_needed = max(target_length - (len_combined + len_others), 40)
        spacers = self._partition_integer(host_needed, 3 if not te4 else 4, min_val=8)

        seq_parts = []
        # Host 1
        seq_parts.append(generate_random_dna(spacers[0], gc_content, self.rng.randint(1, 100000)))

        # TE1 start
        start_te1 = len("".join(seq_parts))
        seq_parts.append(seq_te1[:-overlap_len])
        
        # Overlap region
        start_overlap = len("".join(seq_parts))
        seq_parts.append(overlap_seq)
        end_te1 = len("".join(seq_parts)) - 1
        start_te2 = start_overlap

        # TE2 remainder
        seq_parts.append(seq_te2[overlap_len:])
        end_te2 = len("".join(seq_parts)) - 1

        # Host 2
        seq_parts.append(generate_random_dna(spacers[1], gc_content, self.rng.randint(1, 100000)))

        # TE3
        start_te3 = len("".join(seq_parts))
        seq_parts.append(seq_te3)
        end_te3 = len("".join(seq_parts)) - 1

        # Host 3
        seq_parts.append(generate_random_dna(spacers[2], gc_content, self.rng.randint(1, 100000)))

        if te4:
            start_te4 = len("".join(seq_parts))
            seq_parts.append(seq_te4)
            end_te4 = len("".join(seq_parts)) - 1
            seq_parts.append(generate_random_dna(spacers[3], gc_content, self.rng.randint(1, 100000)))

        final_seq = "".join(seq_parts)

        elements = {}
        elements["TE1_Mariner"] = GroundTruthElement(
            element_id="TE1_Mariner",
            family_name="TE1_Mariner",
            full_start=start_te1,
            full_end=end_te1,
            is_interrupted=False,
            interrupted_by=["TE2_PiggyBac (Overlapped)"],
            nesting_level=0,
            fragments=[GroundTruthFragment("TE1_Mariner", start_te1, end_te1, 0, len(seq_te1) - 1, "overlapping", seq_te1)],
            canonical_sequence=seq_te1
        )
        elements["TE2_PiggyBac"] = GroundTruthElement(
            element_id="TE2_PiggyBac",
            family_name="TE2_PiggyBac",
            full_start=start_te2,
            full_end=end_te2,
            is_interrupted=False,
            interrupted_by=["TE1_Mariner (Overlapped)"],
            nesting_level=0,
            fragments=[GroundTruthFragment("TE2_PiggyBac", start_te2, end_te2, 0, len(seq_te2) - 1, "overlapping", seq_te2)],
            canonical_sequence=seq_te2
        )
        elements["TE3_Gypsy"] = GroundTruthElement(
            element_id="TE3_Gypsy",
            family_name="TE3_Gypsy",
            full_start=start_te3,
            full_end=end_te3,
            is_interrupted=False,
            nesting_level=0,
            fragments=[GroundTruthFragment("TE3_Gypsy", start_te3, end_te3, 0, len(seq_te3) - 1, "intact", seq_te3)],
            canonical_sequence=seq_te3
        )
        if te4:
            elements["TE4_Helitron"] = GroundTruthElement(
                element_id="TE4_Helitron",
                family_name="TE4_Helitron",
                full_start=start_te4,
                full_end=end_te4,
                is_interrupted=False,
                nesting_level=0,
                fragments=[GroundTruthFragment("TE4_Helitron", start_te4, end_te4, 0, len(seq_te4) - 1, "intact", seq_te4)],
                canonical_sequence=seq_te4
            )

        history.append(f"1. TE1_Mariner and TE2_PiggyBac simulated with a {overlap_len} bp coordinate intersection/overlap.")
        history.append("2. TE3_Gypsy placed downstream.")
        if te4:
            history.append("3. TE4_Helitron placed downstream.")

        return SimulationResult(
            dna_sequence=final_seq,
            total_length=len(final_seq),
            scenario_name="Overlapping Transposons (Coordinate Collision)",
            scenario_description=(
                f"Simulates TE1_Mariner and TE2_PiggyBac sharing a {overlap_len}-bp overlapping region "
                f"between coordinates {start_overlap} and {end_te1}."
            ),
            elements=elements,
            reference_library={k: m.full_sequence for k, m in models.items()},
            insertion_history=history
        )

    def _simulate_independent(self,
                             target_length: int,
                             models: Dict[str, TransposonModel],
                             gc_content: float) -> SimulationResult:
        """Simulate independent, non-intersecting transposons (baseline control)."""
        history = []
        selected_keys = list(models.keys())
        te_seqs = [models[k].full_sequence for k in selected_keys]
        total_te_len = sum(len(s) for s in te_seqs)

        host_needed = max(target_length - total_te_len, 30)
        spacers = self._partition_integer(host_needed, len(selected_keys) + 1, min_val=6)

        seq_parts = []
        elements = {}

        for idx, k in enumerate(selected_keys):
            # Spacer
            seq_parts.append(generate_random_dna(spacers[idx], gc_content, self.rng.randint(1, 100000)))
            # TE
            start_te = len("".join(seq_parts))
            seq = models[k].full_sequence
            seq_parts.append(seq)
            end_te = len("".join(seq_parts)) - 1

            elements[k] = GroundTruthElement(
                element_id=k,
                family_name=k,
                full_start=start_te,
                full_end=end_te,
                is_interrupted=False,
                interrupted_by=[],
                nesting_level=0,
                fragments=[GroundTruthFragment(k, start_te, end_te, 0, len(seq) - 1, "intact", seq)],
                canonical_sequence=seq
            )
            history.append(f"{idx+1}. Inserted {k} at continuous span [{start_te}, {end_te}].")

        # Final trailing spacer
        seq_parts.append(generate_random_dna(spacers[-1], gc_content, self.rng.randint(1, 100000)))
        final_seq = "".join(seq_parts)

        return SimulationResult(
            dna_sequence=final_seq,
            total_length=len(final_seq),
            scenario_name="Independent Transposons (Baseline Control)",
            scenario_description="Simulates 3 or 4 transposable elements placed independently with no intersections.",
            elements=elements,
            reference_library={k: m.full_sequence for k, m in models.items()},
            insertion_history=history
        )

    def _partition_integer(self, total: int, parts: int, min_val: int = 5) -> List[int]:
        """Randomly partition an integer into `parts` positive numbers >= min_val."""
        if total < parts * min_val:
            return [max(min_val, total // parts)] * parts
        
        rem = total - parts * min_val
        cuts = sorted([self.rng.randint(0, rem) for _ in range(parts - 1)])
        cuts = [0] + cuts + [rem]
        result = [min_val + (cuts[i+1] - cuts[i]) for i in range(parts)]
        return result
