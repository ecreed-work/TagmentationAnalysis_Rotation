from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, List, Literal, Dict
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from Bio.Align import PairwiseAligner
from collections import Counter
from Bio import SeqIO


@dataclass(frozen=True)
class SpacerHit:
    strand: Literal["+","-"]
    pam_start: int
    pam_end: int
    spacer_start: int
    spacer_end: int
    pam_dna: str
    spacer_top: str   # DNA by default


IUPAC_DNA: Dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]",
}


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]

def load_genome(filename):
    from Bio import SeqIO
    import re
    genome_sequence = str(next(SeqIO.parse(filename, "fasta")).seq)
    genome_sequence = re.sub(r'^>.*\n?','',genome_sequence,flags=re.MULTILINE)
    genome_sequence = genome_sequence.replace('\n','').replace('\r','')
    return genome_sequence


def map_guide_to_genome(guide, genome):
    """
    Map a guide RNA (DNA alphabet) to a genome sequence.

    Parameters
    ----------
    guide : str
        Guide RNA sequence (DNA alphabet: A/C/G/T)
    genome : str
        Genome sequence (DNA alphabet)

    Returns
    -------
    hits : list of dict
        Each hit contains:
        - start: 0-based start index in genome
        - end:   end index (exclusive)
        - strand: '+' or '-'
        - sequence: matched genome sequence
    """
    guide = guide.upper()
    genome = genome.upper()

    rc_guide = reverse_complement(guide)
    guide_len = len(guide)

    hits = []

    # Search forward strand
    start = 0
    while True:
        idx = genome.find(guide, start)
        if idx == -1:
            break
        hits.append({
            "start": idx,
            "end": idx + guide_len,
            "strand": "+",
            "sequence": genome[idx:idx + guide_len]
        })
        start = idx + 1

    # Search reverse strand
    start = 0
    while True:
        idx = genome.find(rc_guide, start)
        if idx == -1:
            break
        hits.append({
            "start": idx,
            "end": idx + guide_len,
            "strand": "-",
            "sequence": genome[idx:idx + guide_len]
        })
        start = idx + 1

    return hits


def extract_sequence_windows(positions, genome, window_bp):
    """
    Extract genomic sequence windows around given positions.

    Parameters
    ----------
    positions : list of int
        Genomic coordinates (0-based).
    genome : str
        Genome sequence.
    window_bp : int
        Number of base pairs to include on each side of the position.

    Returns
    -------
    list of dict
        Each entry contains:
        - position: original genomic position
        - start: window start (clipped)
        - end: window end (exclusive)
        - sequence: extracted sequence
    """
    genome = genome.upper()
    genome_len = len(genome)

    windows = []

    for pos in positions:
        if pos < 0 or pos >= genome_len:
            continue  # skip invalid positions

        start = max(0, pos - window_bp)
        end = min(genome_len, pos + window_bp + 1)

        windows.append({
            "position": pos,
            "start": start,
            "end": end,
            "sequence": genome[start:end]
        })

    return windows

def at_richness(sequence: str, *, ignore_ambiguous: bool = True) -> float:
        """
        Compute AT richness (fraction of A/T bases) for a nucleic-acid sequence.
        
        Parameters
        ----------
        sequence : str
        DNA/RNA sequence. 'U' is treated as 'T'. Case-insensitive.
        ignore_ambiguous : bool, default True
        If True, compute A/T fraction over only unambiguous bases (A/C/G/T).
        If False, compute A/T fraction over all characters except whitespace/gaps.
        
        Returns
        -------
        float
        AT richness in [0, 1]. Returns 0.0 if no valid bases are found.
        """
        if sequence is None:
            raise ValueError("sequence must be a string, not None")
            
        seq = sequence.upper().replace("U", "T")
        
        # Remove common whitespace and gap characters
        
        seq = "".join(ch for ch in seq if ch not in {" ", "\n", "\r", "\t", "-"})
        
        if ignore_ambiguous:
            valid = [ch for ch in seq if ch in {"A", "C", "G", "T"}]
            denom = len(valid)
            at = sum(ch in {"A", "T"} for ch in valid)
        else:
            denom = len(seq)
            at = sum(ch in {"A", "T"} for ch in seq)
            
        return (at / denom) if denom else 0.0

def sliding_window(s, window):
    return [s[i:i+window] for i in range(0, len(s) - window + 1)]

def at_sliding_window( seq, window ):
    import TagmentationAnalysis.postprocessHelpers as hp
    #for window, compute at-richness and store in a list
    chunks = sliding_window( seq, window )
    at = []
    for ii in range(0,len(chunks)):
        at.append( at_richness( chunks[ ii ] ))
    return at

def local_align(query, target):
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -0.5

    aln = aligner.align(query, target)[0]  # best
    return aln

def load_fasta_as_dict(fasta_path: str) -> dict[str, str]:
    seqs = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seqs[rec.id] = str(rec.seq).upper()
    return seqs


def _pam_matches_at(seq: str, i: int, pam: str) -> bool:
    if i < 0 or i + len(pam) > len(seq):
        return False
    for j, code in enumerate(pam.upper()):
        base = seq[i + j].upper().replace("U", "T")
        allowed = IUPAC_DNA[code]
        if allowed.startswith("["):
            if base not in allowed.strip("[]"):
                return False
        else:
            if base != allowed:
                return False
    return True


def extract_spacers(
    sequence: str,
    pam: str,
    spacer_len: int = 20,
    spacer_side: Literal["5prime", "3prime"] = "5prime",
    search_strands: Literal["+", "-", "both"] = "both",
    return_rna: bool = False,   # <-- DNA is now default
) -> List[SpacerHit]:
    """
    Extract all possible spacer sequences adjacent to a PAM.
    by default sequences always returned as top strand sequence, 
    so reverse complement needs to be applied to get actual sequence on bottom strand. 
    Returns DNA spacers by default. Set return_rna=True to return RNA (T->U).
    """
    seq = sequence.replace(" ", "").replace("\n", "")
    pam = pam.upper()

    def dna_or_rna(s: str) -> str:
        return s.replace("T", "U") if return_rna else s

    hits: List[SpacerHit] = []

    # + strand
    def scan_plus():
        for i in range(len(seq) - len(pam) + 1):
            if not _pam_matches_at(seq, i, pam):
                continue

            pam_start, pam_end = i, i + len(pam)
            if spacer_side == "5prime":
                spacer_start, spacer_end = pam_start - spacer_len, pam_start
            else:
                spacer_start, spacer_end = pam_end, pam_end + spacer_len

            if spacer_start < 0 or spacer_end > len(seq):
                continue

            spacer = seq[spacer_start:spacer_end].upper().replace("U", "T")
            pam_seq = seq[pam_start:pam_end].upper().replace("U", "T")

            hits.append(
                SpacerHit(
                    strand="+",
                    pam_start=pam_start,
                    pam_end=pam_end,
                    spacer_start=spacer_start,
                    spacer_end=spacer_end,
                    pam_dna=pam_seq,
                    spacer_top=dna_or_rna(spacer),
                )
            )

    # - strand
    def scan_minus():
        rc = reverse_complement(seq.replace("U", "T"))
        L = len(seq)

        for i in range(len(rc) - len(pam) + 1):
            if not _pam_matches_at(rc, i, pam):
                continue

            pam_start = L - (i + len(pam))
            pam_end = L - i

            if spacer_side == "5prime":
                spacer_start, spacer_end = pam_end, pam_end + spacer_len
            else:
                spacer_start, spacer_end = pam_start - spacer_len, pam_start

            if spacer_start < 0 or spacer_end > L:
                continue

            spacer_plus = seq[spacer_start:spacer_end].upper().replace("U", "T")
            spacer_guide = reverse_complement(spacer_plus)
            pam_seq = seq[pam_start:pam_end].upper().replace("U", "T")

            hits.append(
                SpacerHit(
                    strand="-",
                    pam_start=pam_start,
                    pam_end=pam_end,
                    spacer_start=spacer_start,
                    spacer_end=spacer_end,
                    pam_dna=pam_seq,
                    spacer_top=spacer_plus,
                )
            )

    if search_strands in ("+", "both"):
        scan_plus()
    if search_strands in ("-", "both"):
        scan_minus()

    return hits


def build_log_odds_pwm(counts_df, background=None, pseudocount=1.0, log_base=2):
    """
    Build a log-odds PWM from a counts matrix (positions x alphabet).

    Parameters
    ----------
    counts_df : pd.DataFrame
        Rows = positions (0..L-1), columns = letters (e.g., A,C,G,T), values = counts.
    background : dict or pd.Series or None
        Background probabilities for each letter. If None, uses uniform over columns.
        Example: {"A":0.25,"C":0.25,"G":0.25,"T":0.25}
    pseudocount : float, default 1.0
        Added to each count to avoid zeros (Laplace smoothing).
    log_base : int or float, default 2
        Base for the logarithm. log2 is standard for “bits”.

    Returns
    -------
    pwm_df : pd.DataFrame
        Log-odds PWM with same shape as counts_df:
        pwm_df.loc[i, base] = log_{log_base}( p_i(base) / bg(base) )
    prob_df : pd.DataFrame
        Smoothed position-specific probabilities used to compute the PWM.
    """
    alphabet = list(counts_df.columns)

    # Background
    if background is None:
        bg = pd.Series({a: 1.0 / len(alphabet) for a in alphabet}, dtype=float)
    else:
        bg = pd.Series(background, dtype=float).reindex(alphabet)
        if bg.isna().any():
            missing = bg[bg.isna()].index.tolist()
            raise ValueError(f"Background missing letters: {missing}")
        if (bg <= 0).any():
            raise ValueError("Background probabilities must be > 0")

    # Smooth counts -> probabilities
    smoothed = counts_df.astype(float) + float(pseudocount)
    prob_df = smoothed.div(smoothed.sum(axis=1), axis=0)

    # Log-odds
    denom = bg.values  # length K
    ratio = prob_df.values / denom  # (L x K)

    if log_base == 2:
        pwm = np.log2(ratio)
    else:
        pwm = np.log(ratio) / np.log(log_base)

    pwm_df = pd.DataFrame(pwm, index=counts_df.index, columns=alphabet)
    return pwm_df, prob_df

def score_sequence_pwm(seq, pwm_df, *, unknown_penalty=None):
    """
    Score a sequence with a log-odds PWM (sum across positions).

    Parameters
    ----------
    seq : str
        Sequence to score (length must equal #rows in pwm_df).
    pwm_df : pd.DataFrame
        Log-odds PWM (positions x alphabet).
    unknown_penalty : float or None, default None
        What to do if a character is not in pwm_df.columns:
        - None: raise ValueError
        - float: add this value for unknown characters

    Returns
    -------
    float
        Total log-odds score.
    """
    seq = seq.upper()
    L = pwm_df.shape[0]
    if len(seq) != L:
        raise ValueError(f"Sequence length {len(seq)} != PWM length {L}")

    col_index = {c: j for j, c in enumerate(pwm_df.columns)}
    vals = pwm_df.values

    total = 0.0
    for i, ch in enumerate(seq):
        j = col_index.get(ch)
        if j is None:
            if unknown_penalty is None:
                raise ValueError(f"Unknown base '{ch}' at position {i}")
            total += float(unknown_penalty)
        else:
            total += float(vals[i, j])

    return total

import numpy as np

def pwm_max_possible_score(pwm_df):
    """
    Maximum possible score achievable by this PWM over any sequence of length L,
    i.e., sum of the per-position best letter log-odds.
    """
    return float(pwm_df.max(axis=1).sum())

def predict_insertion_from_spacer(guideRNA, genome_sequence):
    #from spacer find insertions that occur at the 'correct' distance: 
    distance = 71
    window = 10
    
    grna = map_guide_to_genome(guideRNA,genome_sequence)
    if grna[0]['strand'] == '+':
        predicted_pos = grna[0]['end']+3+distance
    else:
        predicted_pos = grna[0]['start']-3-distance

    df = pd.DataFrame(columns=["PAM","protospacer","gRNA_strand","intervening","insertion_length","target_site",
                               "PAM_start","PAM_end","proto_start","proto_end","aln_score",
                               "intervene_start","intervene_end","ts_start","ts_end"])

    df['proto_start']= [ grna[0]['start'] ]
    df['proto_end']= [ grna[0]['end'] ]
    df['protospacer']= [ grna[0]['sequence'] ]
    if grna[0]['strand'] == '+':
        df['PAM_start'] = [ grna[0]['end'] ]
        df['PAM_end']= [ grna[0]['end']+3 ]
        df['PAM']= [ genome_sequence[ grna[0]['end'] : grna[0]['end']+3 ] ]
        df['intervene_start']= [ grna[0]['end']+3 ]
        df['intervene_end']= [ predicted_pos ]
        df['intervening']= [ genome_sequence[ grna[0]['end']+3 :  predicted_pos  ] ]
        df['target_site']= [ genome_sequence[predicted_pos-window:predicted_pos+window] ]
        df['insertion_length']= [ len(genome_sequence[ grna[0]['end']+3 :  predicted_pos  ]) ]

    else:
        df['PAM_start'] = [ grna[0]['start']-3 ]
        df['PAM_end']= [ grna[0]['start'] ]
        df['PAM']= [ genome_sequence[grna[0]['start']-3:grna[0]['start']]  ] #always list top strand
        df['intervene_start']= [ predicted_pos ]
        df['intervene_end']= [ grna[0]['start']-3 ]
        df['intervening']= [ genome_sequence[predicted_pos:grna[0]['start']-3] ]
        df['target_site']= [ genome_sequence[predicted_pos-window:predicted_pos+window]  ]
        df['insertion_length']=[ len(genome_sequence[predicted_pos:grna[0]['start']-3] ) ]

    df['ts_start']= [ predicted_pos-window ]
    df['ts_end']= [ predicted_pos+window ]
    df['gRNA_strand']= [ grna[0]['strand'] ]
    return df

def degenerate_hamming(window: str, motif_iupac: str) -> int:
    """
    Returns number of mismatches under IUPAC degeneracy.
    0 = perfect match, higher = worse.
    """
    window = window.upper().replace("U", "T")
    motif_iupac = motif_iupac.upper()
    if len(window) != len(motif_iupac):
        raise ValueError("window and motif must have the same length")

    mismatches = 0
    for b, code in zip(window, motif_iupac):
        allowed = IUPAC_DNA.get(code)
        if allowed is None:
            raise ValueError(f"Unknown IUPAC code: {code}")
        if b not in allowed:
            mismatches += 1
    return mismatches

def degenerate_match_count(window: str, motif_iupac: str) -> int:
    """Returns number of matching positions (higher = better)."""
    return len(motif_iupac) - degenerate_hamming(window, motif_iupac)

def scan_degenerate_hamming(seq: str, motif_iupac: str):
    seq = seq.upper().replace("U", "T")
    m = len(motif_iupac)
    hits = []
    for i in range(len(seq) - m + 1):
        w = seq[i:i+m]
        d = degenerate_hamming(w, motif_iupac)
        hits.append((i, d, w))
    # sort by best (fewest mismatches)
    hits.sort(key=lambda t: t[1])
    return hits

def best_hit_fewest_mismatches_closest_center(seq: str, motif_iupac: str):
    """
    Returns the (start, mismatches, window) hit with:
      1) minimum mismatches
      2) closest window-center to sequence center
      3) (tie-break) smallest start index
    """
    seq = seq.upper().replace("U", "T")
    m = len(motif_iupac)

    seq_center = (len(seq) - 1) / 2.0  # center coordinate of the full sequence

    best = None
    best_key = None

    for start in range(len(seq) - m + 1):
        window = seq[start:start + m]
        mismatches = degenerate_hamming(window, motif_iupac)

        window_center = start + (m - 1) / 2.0
        dist_to_center = abs(window_center - seq_center)

        key = (mismatches, dist_to_center, start)  # lexicographic sort: smaller is better

        if best_key is None or key < best_key:
            best_key = key
            best = (start, mismatches, window)

    return best

def _iupac_match(seq_window: str, pattern: str) -> bool:
    """Return True if seq_window matches IUPAC pattern (same length)."""
    if len(seq_window) != len(pattern):
        return False
    seq_window = seq_window.upper().replace("U", "T")
    pattern = pattern.upper().replace("U", "T")
    for b, code in zip(seq_window, pattern):
        allowed = IUPAC_DNA.get(code)
        if allowed is None:
            raise ValueError(f"Unknown IUPAC code in PAM motif: {code}")
        if b not in allowed:
            return False
    return True


def find_pam_sites_with_optional_protospacer(
    segment: str,
    pam_motif: str,
    protospacer: Optional[str] = None,
    pam_side: str = "3prime",                 # "3prime"/"3'" or "5prime"/"5'"
    include_reverse_complement: bool = False, # scan both strands
) -> List[Dict]:
    """
    Find PAM sites in `segment` using an IUPAC-degenerate PAM motif.
    If `protospacer` is provided, require an adjacent exact protospacer match.

    pam_side interpretation (on the scanned strand):
      - "3prime": protospacer immediately upstream of PAM (protospacer + PAM)
      - "5prime": protospacer immediately downstream of PAM (PAM + protospacer)

    Returns a list of dict hits with coordinates in the + orientation of `segment`.
    """
    # normalize
    segment = segment.upper().replace("U", "T")
    pam_motif = pam_motif.upper().replace("U", "T")
    if protospacer is not None:
        protospacer = protospacer.upper().replace("U", "T")

    # normalize pam_side
    ps = pam_side.strip().lower()
    if ps in {"3", "3prime", "3'"}:
        ps = "3prime"
    elif ps in {"5", "5prime", "5'"}:
        ps = "5prime"
    else:
        raise ValueError("pam_side must be '3prime'/'3'' or '5prime'/'5'' (or '3'/'5').")

    pam_len = len(pam_motif)
    if pam_len == 0:
        raise ValueError("pam_motif cannot be empty.")

    proto_len = len(protospacer) if protospacer is not None else 0

    hits: List[Dict] = []

    def scan(seq: str, strand: str):
        # seq is the sequence in the orientation being scanned
        for pam_start in range(0, len(seq) - pam_len + 1):
            pam_window = seq[pam_start:pam_start + pam_len]
            if not _iupac_match(pam_window, pam_motif):
                continue

            # If no protospacer requirement, record PAM hit
            if protospacer is None:
                # Map to + coords
                pos_plus = pam_start if strand == "+" else (len(segment) - pam_len - pam_start)
                hits.append({
                    "pam_start": pos_plus,
                    "pam_end": pos_plus + pam_len,   # end-exclusive in + coords
                    "strand": strand,
                    "pam_seq": pam_window,
                })
                continue

            # Determine where protospacer must lie relative to PAM
            if ps == "3prime":
                # protospacer + PAM => protospacer immediately upstream
                proto_start = pam_start - proto_len
                proto_end = pam_start
            else:
                # PAM + protospacer => protospacer immediately downstream
                proto_start = pam_start + pam_len
                proto_end = proto_start + proto_len

            if proto_start < 0 or proto_end > len(seq):
                continue

            proto_window = seq[proto_start:proto_end]
            if proto_window != protospacer:
                continue

            # Map both PAM and protospacer coordinates to + orientation of original segment
            if strand == "+":
                pam_start_plus = pam_start
                proto_start_plus = proto_start
            else:
                # For RC scanning, a window starting at i maps to + start = L - window_len - i
                L = len(segment)
                pam_start_plus = L - pam_len - pam_start
                proto_start_plus = L - proto_len - proto_start

            hits.append({
                "pam_start": pam_start_plus,
                "pam_end": pam_start_plus + pam_len,              # end-exclusive
                "protospacer_start": proto_start_plus,
                "protospacer_end": proto_start_plus + proto_len,  # end-exclusive
                "strand": strand,
                "pam_seq": pam_window,
                "protospacer_seq": proto_window,
                "pam_side": ps,
            })

    # scan + strand
    scan(segment, "+")

    # scan - strand (optional)
    if include_reverse_complement:
        scan(reverse_complement(segment), "-")

    # sort by PAM start then strand
    hits.sort(key=lambda d: (d["pam_start"], d["strand"]))
    return hits
