"""
detector.py - Transposable Element Detection Engine with Intersection Resolution
================================================================================
Implements multiple detection algorithms to identify transposable element
positions (start, end) and analyzes how algorithms behave when transposons intersect
(nested insertions, fragmented elements, and coordinate overlaps).

Implemented Algorithms:
1. Split-Aware Synteny & Nesting Matcher (Advanced Chaining Resolver)
2. Naive Continuous Local Alignment Matcher
3. De Novo Structural Matcher (TIR & TSD Scanner)
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Any
import math


def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
            'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))


def smith_waterman_local_align(query: str, ref: str,
                               match_score: int = 2,
                               mismatch_penalty: int = -2,
                               gap_penalty: int = -3) -> Tuple[int, int, int, int, int, str, str]:
    """
    Standard Smith-Waterman local alignment between query and reference.
    Returns:
    (max_score, q_start, q_end, ref_start, ref_end, aligned_q, aligned_ref)
    """
    m, n = len(query), len(ref)
    if m == 0 or n == 0:
        return 0, 0, 0, 0, 0, "", ""

    # Score matrix
    H = [[0] * (n + 1) for _ in range(m + 1)]
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m + 1):
        q_char = query[i - 1]
        for j in range(1, n + 1):
            r_char = ref[j - 1]
            diag = H[i - 1][j - 1] + (match_score if q_char == r_char else mismatch_penalty)
            delete = H[i - 1][j] + gap_penalty
            insert = H[i][j - 1] + gap_penalty
            score = max(0, diag, delete, insert)
            H[i][j] = score
            if score > max_score:
                max_score = score
                max_pos = (i, j)

    if max_score == 0:
        return 0, 0, 0, 0, 0, "", ""

    # Traceback
    i, j = max_pos
    q_end = i - 1
    ref_end = j - 1

    aligned_q = []
    aligned_ref = []

    while i > 0 and j > 0 and H[i][j] > 0:
        curr = H[i][j]
        q_char = query[i - 1]
        r_char = ref[j - 1]
        diag = H[i - 1][j - 1] + (match_score if q_char == r_char else mismatch_penalty)
        delete = H[i - 1][j] + gap_penalty
        insert = H[i][j - 1] + gap_penalty

        if curr == diag:
            aligned_q.append(q_char)
            aligned_ref.append(r_char)
            i -= 1
            j -= 1
        elif curr == delete:
            aligned_q.append(q_char)
            aligned_ref.append('-')
            i -= 1
        elif curr == insert:
            aligned_q.append('-')
            aligned_ref.append(r_char)
            j -= 1
        else:
            break

    q_start = i
    ref_start = j
    aligned_q_str = "".join(reversed(aligned_q))
    aligned_ref_str = "".join(reversed(aligned_ref))

    return max_score, q_start, q_end, ref_start, ref_end, aligned_q_str, aligned_ref_str


def find_all_local_hits(query: str, ref: str,
                        min_match_len: int = 15,
                        kmer_size: int = 10) -> List[Dict[str, Any]]:
    """
    Find all high-scoring local matching fragments of a reference element within the query.
    Uses k-mer seeding followed by extension.
    """
    hits = []
    query_len = len(query)
    ref_len = len(ref)

    # Index k-mers in query
    query_kmers: Dict[str, List[int]] = {}
    for i in range(query_len - kmer_size + 1):
        kmer = query[i:i + kmer_size]
        query_kmers.setdefault(kmer, []).append(i)

    # Scan reference against query k-mers
    raw_seeds = []
    for j in range(ref_len - kmer_size + 1):
        kmer = ref[j:j + kmer_size]
        if kmer in query_kmers:
            for q_idx in query_kmers[kmer]:
                raw_seeds.append((q_idx, j))

    if not raw_seeds:
        # Fall back to Smith-Waterman
        score, qs, qe, rs, re, aq, ar = smith_waterman_local_align(query, ref)
        if (qe - qs + 1) >= min_match_len:
            hits.append({
                "query_start": qs,
                "query_end": qe,
                "ref_start": rs,
                "ref_end": re,
                "length": qe - qs + 1,
                "score": score,
                "aligned_query": aq,
                "aligned_ref": ar
            })
        return hits

    # Group seeds into diagonal runs (ungapped extension)
    diagonal_runs: Dict[int, List[Tuple[int, int]]] = {}
    for q_idx, r_idx in raw_seeds:
        diag = q_idx - r_idx
        diagonal_runs.setdefault(diag, []).append((q_idx, r_idx))

    merged_hits = []
    for diag, seed_list in diagonal_runs.items():
        seed_list.sort(key=lambda x: x[0])
        # Extend runs
        cur_q_start, cur_r_start = seed_list[0]
        cur_q_end = cur_q_start + kmer_size - 1
        cur_r_end = cur_r_start + kmer_size - 1

        for q_pos, r_pos in seed_list[1:]:
            if q_pos <= cur_q_end + 1:
                cur_q_end = max(cur_q_end, q_pos + kmer_size - 1)
                cur_r_end = max(cur_r_end, r_pos + kmer_size - 1)
            else:
                if (cur_q_end - cur_q_start + 1) >= min_match_len:
                    merged_hits.append((cur_q_start, cur_q_end, cur_r_start, cur_r_end))
                cur_q_start, cur_r_start = q_pos, r_pos
                cur_q_end = cur_q_start + kmer_size - 1
                cur_r_end = cur_r_start + kmer_size - 1

        if (cur_q_end - cur_q_start + 1) >= min_match_len:
            merged_hits.append((cur_q_start, cur_q_end, cur_r_start, cur_r_end))

    # Further expand hits left and right
    for qs, qe, rs, re in merged_hits:
        # Extend left
        while qs > 0 and rs > 0 and query[qs - 1] == ref[rs - 1]:
            qs -= 1
            rs -= 1
        # Extend right
        while qe < query_len - 1 and re < ref_len - 1 and query[qe + 1] == ref[re + 1]:
            qe += 1
            re += 1

        match_sub_q = query[qs:qe + 1]
        match_sub_r = ref[rs:re + 1]
        score = sum(2 if a == b else -2 for a, b in zip(match_sub_q, match_sub_r))

        hits.append({
            "query_start": qs,
            "query_end": qe,
            "ref_start": rs,
            "ref_end": re,
            "length": qe - qs + 1,
            "score": score,
            "aligned_query": match_sub_q,
            "aligned_ref": match_sub_r
        })

    # Deduplicate overlapping hits for the same reference
    hits.sort(key=lambda h: (h["query_start"], -h["length"]))
    unique_hits = []
    for h in hits:
        overlap = False
        for u in unique_hits:
            if abs(h["query_start"] - u["query_start"]) <= 2 and abs(h["query_end"] - u["query_end"]) <= 2:
                overlap = True
                break
        if not overlap:
            unique_hits.append(h)

    return unique_hits


@dataclass
class DetectedFragment:
    """A distinct continuous sequence match belonging to a detected element."""
    start: int                  # 0-based inclusive
    end: int                    # 0-based inclusive
    ref_start: int
    ref_end: int
    length: int
    identity: float
    fragment_type: str          # e.g., '5_prime_flank', '3_prime_flank', 'full_intact'
    subsequence: str


@dataclass
class DetectedElement:
    """A predicted transposable element in the artificial genome."""
    element_id: str
    family_name: str
    start: int                  # 0-based inclusive overall outer start
    end: int                    # 0-based inclusive overall outer end
    length: int                 # Outer span length
    is_interrupted: bool = False
    interrupted_by: List[str] = field(default_factory=list)
    parent_element: Optional[str] = None
    nesting_level: int = 0
    confidence_score: float = 1.0
    detection_method: str = ""  # 'split_aware', 'naive_continuous', 'structural_tir'
    fragments: List[DetectedFragment] = field(default_factory=list)
    tsd_motif: Optional[str] = None
    tsd_left_pos: Optional[int] = None
    tsd_right_pos: Optional[int] = None
    behavior_notes: str = ""    # Educational commentary on algorithm behavior


@dataclass
class BenchmarkComparison:
    """Comparison metrics between detected element and ground truth."""
    element_id: str
    ground_truth_start: int
    ground_truth_end: int
    predicted_start: Optional[int]
    predicted_end: Optional[int]
    delta_start: Optional[int]
    delta_end: Optional[int]
    iou_score: float            # Intersection over Union (Jaccard)
    is_detected: bool
    status: str                 # 'Exact Match', 'Split Resolved', 'Truncated/Partial', 'Missed', 'Boundary Shift'
    gt_interrupted: bool
    pred_interrupted: bool
    algorithm_behavior_comment: str


class SplitAwareTransposonDetector:
    """
    Advanced Split-Aware Synteny & Nesting Matcher.
    Capable of identifying intact transposable elements as well as interrupted/nested elements
    where one transposon has jumped inside another, fracturing the parent transposon into fragments.
    """

    def __init__(self, te_library: Dict[str, str], min_fragment_len: int = 14):
        self.te_library = te_library  # Dict of {name: canonical_seq}
        self.min_fragment_len = min_fragment_len

    def detect(self, dna_sequence: str) -> List[DetectedElement]:
        """
        Detect transposable elements in the sequence, reconstructing interrupted
        elements and building the nesting hierarchy.
        """
        # Step 1: Collect all high-scoring fragment alignments for each reference TE
        raw_element_hits: Dict[str, List[Dict[str, Any]]] = {}
        for te_name, ref_seq in self.te_library.items():
            hits = find_all_local_hits(dna_sequence, ref_seq, min_match_len=self.min_fragment_len)
            if hits:
                raw_element_hits[te_name] = sorted(hits, key=lambda x: x["query_start"])

        # Step 2: For each element, evaluate if it is continuous or split/interrupted
        detected_elements: List[DetectedElement] = []

        # First pass: identify all candidate elements and their fragments
        candidate_map: Dict[str, DetectedElement] = {}

        for te_name, hits in raw_element_hits.items():
            ref_seq = self.te_library[te_name]
            ref_len = len(ref_seq)

            if len(hits) == 1:
                # Single continuous hit
                h = hits[0]
                frag = DetectedFragment(
                    start=h["query_start"],
                    end=h["query_end"],
                    ref_start=h["ref_start"],
                    ref_end=h["ref_end"],
                    length=h["length"],
                    identity=1.0,
                    fragment_type="full_intact",
                    subsequence=dna_sequence[h["query_start"]:h["query_end"] + 1]
                )
                cand = DetectedElement(
                    element_id=te_name,
                    family_name=te_name,
                    start=h["query_start"],
                    end=h["query_end"],
                    length=h["length"],
                    is_interrupted=False,
                    interrupted_by=[],
                    parent_element=None,
                    nesting_level=0,
                    confidence_score=round(h["length"] / ref_len, 3),
                    detection_method="Split-Aware Matcher",
                    fragments=[frag],
                    behavior_notes="Contiguous unbroken alignment across entire element."
                )
                candidate_map[te_name] = cand

            elif len(hits) >= 2:
                # Multiple hits for the same reference element: Check for split/interruption synteny
                # Colinearity check: 5' piece and 3' piece in reference order
                h_first = hits[0]
                h_last = hits[-1]
                
                # Check if fragments are colinear along the reference element
                is_colinear = (h_first["ref_start"] < h_last["ref_start"] and
                               h_first["query_start"] < h_last["query_start"])
                
                intervening_gap = h_last["query_start"] - (h_first["query_end"] + 1)
                
                if is_colinear and intervening_gap >= 15:
                    # Clear nested/interrupted signature!
                    frags = []
                    total_cov = 0
                    for idx, h in enumerate(hits):
                        ftype = "5'_fragment" if idx == 0 else ("3'_fragment" if idx == len(hits) - 1 else f"internal_frag_{idx}")
                        frags.append(DetectedFragment(
                            start=h["query_start"],
                            end=h["query_end"],
                            ref_start=h["ref_start"],
                            ref_end=h["ref_end"],
                            length=h["length"],
                            identity=1.0,
                            fragment_type=ftype,
                            subsequence=dna_sequence[h["query_start"]:h["query_end"] + 1]
                        ))
                        total_cov += h["length"]

                    outer_start = h_first["query_start"]
                    outer_end = h_last["query_end"]
                    span_len = outer_end - outer_start + 1

                    cand = DetectedElement(
                        element_id=te_name,
                        family_name=te_name,
                        start=outer_start,
                        end=outer_end,
                        length=span_len,
                        is_interrupted=True,
                        interrupted_by=[],  # Will be populated in Step 3
                        parent_element=None,
                        nesting_level=0,
                        confidence_score=round(total_cov / ref_len, 3),
                        detection_method="Split-Aware Matcher",
                        fragments=frags,
                        behavior_notes=(
                            f"Identified as INTERRUPTED/NESTED! Reconstructed across {intervening_gap} bp internal gap. "
                            f"Chained {len(frags)} collinear fragments (5' [{h_first['query_start']}..{h_first['query_end']}] "
                            f"and 3' [{h_last['query_start']}..{h_last['query_end']}])."
                        )
                    )
                    candidate_map[te_name] = cand
                else:
                    # Multiple non-colinear hits; select best hit
                    best_hit = max(hits, key=lambda x: x["length"])
                    frag = DetectedFragment(
                        start=best_hit["query_start"],
                        end=best_hit["query_end"],
                        ref_start=best_hit["ref_start"],
                        ref_end=best_hit["ref_end"],
                        length=best_hit["length"],
                        identity=1.0,
                        fragment_type="partial_best",
                        subsequence=dna_sequence[best_hit["query_start"]:best_hit["query_end"] + 1]
                    )
                    cand = DetectedElement(
                        element_id=te_name,
                        family_name=te_name,
                        start=best_hit["query_start"],
                        end=best_hit["query_end"],
                        length=best_hit["length"],
                        is_interrupted=False,
                        parent_element=None,
                        confidence_score=round(best_hit["length"] / ref_len, 3),
                        detection_method="Split-Aware Matcher",
                        fragments=[frag],
                        behavior_notes="Multiple matches detected, reported largest hit."
                    )
                    candidate_map[te_name] = cand

        # Step 3: Build Nesting / Interruption Relationships
        # If element B is located entirely inside the span of element A [start_A, end_A],
        # then B is nested inside A!
        elements_list = list(candidate_map.values())
        # Sort by outer span descending so outermost parent elements come first
        elements_list.sort(key=lambda e: (e.end - e.start), reverse=True)

        for i, parent in enumerate(elements_list):
            for child in elements_list[i + 1:]:
                # Check if child is strictly inside parent's outer boundaries
                if parent.start < child.start and child.end < parent.end:
                    child.parent_element = parent.element_id
                    child.nesting_level = parent.nesting_level + 1
                    if child.element_id not in parent.interrupted_by:
                        parent.interrupted_by.append(child.element_id)

        # Update notes with nesting hierarchy
        for e in elements_list:
            if e.is_interrupted and e.interrupted_by:
                nested_children = ", ".join(e.interrupted_by)
                e.behavior_notes += f" Encloses nested insertion(s): [{nested_children}]."
            if e.parent_element:
                e.behavior_notes += f" Nested inside parent '{e.parent_element}' (Nesting Level {e.nesting_level})."

        # Step 4: Scan for flanking Target Site Duplications (TSDs)
        for e in elements_list:
            self._detect_flanking_tsds(dna_sequence, e)

        # Sort by genomic start coordinate
        elements_list.sort(key=lambda e: e.start)
        return elements_list

    def _detect_flanking_tsds(self, seq: str, elem: DetectedElement):
        """Search for short target site duplications immediately flanking the element."""
        # Check candidate TSD lengths: 2, 3, 4, 5 bp
        best_tsd = None
        for k in [4, 3, 2, 5]:
            left_tsd_cand = seq[max(0, elem.start - k): elem.start]
            right_tsd_cand = seq[elem.end + 1: elem.end + 1 + k]
            if len(left_tsd_cand) == k and left_tsd_cand == right_tsd_cand:
                best_tsd = (left_tsd_cand, elem.start - k, elem.end + 1)
                break

        if best_tsd:
            elem.tsd_motif = best_tsd[0]
            elem.tsd_left_pos = best_tsd[1]
            elem.tsd_right_pos = best_tsd[2]


class NaiveContinuousTransposonDetector:
    """
    Naive Continuous Local Alignment Matcher.
    Assumes all transposable elements are continuous unbroken segments.
    Does not perform fragment synteny chaining or nesting hierarchy detection.
    Illustrates algorithm failure modes during transposon intersections!
    """

    def __init__(self, te_library: Dict[str, str], min_score_threshold: int = 25):
        self.te_library = te_library
        self.min_score_threshold = min_score_threshold

    def detect(self, dna_sequence: str) -> List[DetectedElement]:
        """
        Detect transposable elements assuming continuous, unbroken alignment.
        """
        detected = []

        for te_name, ref_seq in self.te_library.items():
            ref_len = len(ref_seq)
            score, qs, qe, rs, re, aq, ar = smith_waterman_local_align(dna_sequence, ref_seq)

            if score >= self.min_score_threshold:
                match_len = qe - qs + 1
                coverage = match_len / ref_len

                is_truncated = coverage < 0.85
                notes = ""
                if is_truncated:
                    notes = (
                        f"NAIVE BEHAVIOR WARNING: Detected truncated segment of length {match_len} bp "
                        f"({coverage:.1%} coverage of reference {ref_len} bp). "
                        f"Because this algorithm only seeks continuous alignments, it misses the other half "
                        f"of an interrupted element and fails to recognize nesting!"
                    )
                else:
                    notes = "Detected complete unbroken continuous alignment."

                frag = DetectedFragment(
                    start=qs,
                    end=qe,
                    ref_start=rs,
                    ref_end=re,
                    length=match_len,
                    identity=1.0,
                    fragment_type="naive_continuous",
                    subsequence=dna_sequence[qs:qe + 1]
                )

                elem = DetectedElement(
                    element_id=te_name,
                    family_name=te_name,
                    start=qs,
                    end=qe,
                    length=match_len,
                    is_interrupted=False,       # Naive matcher is unaware of interruptions!
                    interrupted_by=[],
                    parent_element=None,
                    nesting_level=0,
                    confidence_score=round(coverage, 3),
                    detection_method="Naive Continuous Matcher",
                    fragments=[frag],
                    behavior_notes=notes
                )
                detected.append(elem)

        detected.sort(key=lambda e: e.start)
        return detected


class StructuralDenovoTransposonDetector:
    """
    De Novo Structural Transposon Detector.
    Searches for transposable elements purely from structural hallmarks:
    Terminal Inverted Repeats (TIRs) and Target Site Duplications (TSDs).
    Requires NO reference sequence library!
    """

    def __init__(self, min_tir_len: int = 7, max_tir_len: int = 12,
                 min_te_len: int = 35, max_te_len: int = 250):
        self.min_tir_len = min_tir_len
        self.max_tir_len = max_tir_len
        self.min_te_len = min_te_len
        self.max_te_len = max_te_len

    def detect(self, dna_sequence: str) -> List[DetectedElement]:
        """
        Scan for paired terminal inverted repeats (TIRs) with optional flanking TSDs.
        """
        seq_len = len(dna_sequence)
        candidates: List[DetectedElement] = []
        element_counter = 1

        for tir_len in range(self.max_tir_len, self.min_tir_len - 1, -1):
            for i in range(seq_len - self.min_te_len):
                tir_5 = dna_sequence[i:i + tir_len]
                tir_3_expected = reverse_complement(tir_5)

                # Search downstream within allowed TE lengths
                min_j = i + self.min_te_len - tir_len
                max_j = min(i + self.max_te_len - tir_len, seq_len - tir_len)

                for j in range(min_j, max_j + 1):
                    tir_3_actual = dna_sequence[j:j + tir_len]
                    if tir_3_actual == tir_3_expected:
                        te_start = i
                        te_end = j + tir_len - 1
                        total_len = te_end - te_start + 1

                        # Check flanking TSD (2 to 4 bp)
                        tsd_found = None
                        for tsd_k in [4, 3, 2]:
                            if te_start >= tsd_k and te_end + tsd_k < seq_len:
                                left_tsd = dna_sequence[te_start - tsd_k:te_start]
                                right_tsd = dna_sequence[te_end + 1:te_end + 1 + tsd_k]
                                if left_tsd == right_tsd:
                                    tsd_found = left_tsd
                                    break

                        # Avoid duplicate hits overlapping heavily with existing candidates
                        overlap = False
                        for c in candidates:
                            if abs(c.start - te_start) <= 3 and abs(c.end - te_end) <= 3:
                                overlap = True
                                break

                        if not overlap:
                            elem_id = f"Denovo_TIR_TE_{element_counter}"
                            element_counter += 1
                            notes = (
                                f"De Novo Discovery: Paired {tir_len}-bp TIRs (5'-{tir_5}...{tir_3_actual}-3'). "
                                f"TSD motif: {tsd_found or 'None'}."
                            )

                            frag = DetectedFragment(
                                start=te_start,
                                end=te_end,
                                ref_start=0,
                                ref_end=total_len - 1,
                                length=total_len,
                                identity=1.0,
                                fragment_type="denovo_tir_pair",
                                subsequence=dna_sequence[te_start:te_end + 1]
                            )

                            elem = DetectedElement(
                                element_id=elem_id,
                                family_name="Structural_TIR_Family",
                                start=te_start,
                                end=te_end,
                                length=total_len,
                                is_interrupted=False,
                                confidence_score=0.90 if tsd_found else 0.75,
                                detection_method="Structural De Novo (TIR/TSD)",
                                fragments=[frag],
                                tsd_motif=tsd_found,
                                behavior_notes=notes
                            )
                            candidates.append(elem)

        # Detect nesting among structural candidates (outer TIR pair enclosing inner TIR pair)
        candidates.sort(key=lambda e: (e.end - e.start), reverse=True)
        for i, parent in enumerate(candidates):
            for child in candidates[i + 1:]:
                if parent.start < child.start and child.end < parent.end:
                    child.parent_element = parent.element_id
                    child.nesting_level = parent.nesting_level + 1
                    parent.is_interrupted = True
                    if child.element_id not in parent.interrupted_by:
                        parent.interrupted_by.append(child.element_id)
                        parent.behavior_notes += f" Outer TIRs enclose nested internal element {child.element_id}."

        candidates.sort(key=lambda e: e.start)
        return candidates


def evaluate_detection(ground_truth_elements: Dict[str, Any],
                       detected_elements: List[DetectedElement]) -> List[BenchmarkComparison]:
    """
    Evaluate predicted transposable elements against ground truth.
    Calculates Intersection-Over-Union (IoU), boundary errors (delta start/end),
    and analyzes specific behavioral traits on intersected elements.
    """
    evaluations: List[BenchmarkComparison] = []

    # Map predictions by best IoU match to each ground truth element
    for gt_id, gt_elem in ground_truth_elements.items():
        gt_start = gt_elem.full_start
        gt_end = gt_elem.full_end
        gt_len = gt_end - gt_start + 1

        # Find best matching prediction
        best_pred = None
        best_iou = 0.0

        for pred in detected_elements:
            # Check overlap
            inter_start = max(gt_start, pred.start)
            inter_end = min(gt_end, pred.end)
            if inter_end >= inter_start:
                inter_len = inter_end - inter_start + 1
                union_len = (gt_len + pred.length) - inter_len
                iou = inter_len / union_len if union_len > 0 else 0.0
                if iou > best_iou:
                    best_iou = iou
                    best_pred = pred

        if best_pred and best_iou >= 0.20:
            d_start = abs(best_pred.start - gt_start)
            d_end = abs(best_pred.end - gt_end)

            if d_start == 0 and d_end == 0:
                if gt_elem.is_interrupted:
                    status = "Exact Match (Interruption Resolved)"
                else:
                    status = "Exact Match (Continuous)"
            elif d_start <= 3 and d_end <= 3:
                status = "Near Exact Match (Boundary Shift <= 3bp)"
            elif best_iou >= 0.50:
                status = "Partial / Truncated Match"
            else:
                status = "Fragmented / Low Overlap"

            comment = (
                f"Predicted [{best_pred.start}..{best_pred.end}] vs GT [{gt_start}..{gt_end}]. "
                f"IoU: {best_iou:.3f}. Delta Start: {d_start} bp, Delta End: {d_end} bp. "
                f"Method: {best_pred.detection_method}. {best_pred.behavior_notes}"
            )

            evaluations.append(BenchmarkComparison(
                element_id=gt_id,
                ground_truth_start=gt_start,
                ground_truth_end=gt_end,
                predicted_start=best_pred.start,
                predicted_end=best_pred.end,
                delta_start=d_start,
                delta_end=d_end,
                iou_score=round(best_iou, 3),
                is_detected=True,
                status=status,
                gt_interrupted=gt_elem.is_interrupted,
                pred_interrupted=best_pred.is_interrupted,
                algorithm_behavior_comment=comment
            ))
        else:
            evaluations.append(BenchmarkComparison(
                element_id=gt_id,
                ground_truth_start=gt_start,
                ground_truth_end=gt_end,
                predicted_start=None,
                predicted_end=None,
                delta_start=None,
                delta_end=None,
                iou_score=0.0,
                is_detected=False,
                status="Missed / Undetected",
                gt_interrupted=gt_elem.is_interrupted,
                pred_interrupted=False,
                algorithm_behavior_comment="No corresponding predicted element matched this ground truth region."
            ))

    return evaluations
