from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, List, Literal, Dict
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from Bio.Align import PairwiseAligner
from collections import Counter
import TagmentationAnalysis.postprocessHelpers as hp
import statistics

def which(df, pos): #operates on tsv map
    return df[df["ins0"] == pos]

def between(df,s,e): #operates on tsv map
    subset = df[df["ins0"].between(s, e)]
    return subset


def percent_within_distance(insertions, target_position, window_bp):
    """
    Compute the percentage of insertions within ±window_bp of a target position.

    Parameters
    ----------
    insertions : iterable of int
        Genomic coordinates of insertions (e.g., start positions).
    target_position : int
        Genomic coordinate to compare against.
    window_bp : int
        Distance window in base pairs (±window_bp).

    Returns
    -------
    percent : float
        Percentage of insertions within the window.
    """
    if len(insertions) == 0:
        return 0.0

    count_within = sum(
        abs(pos - target_position) <= window_bp
        for pos in insertions
    )

    return 100.0 * count_within / len(insertions)

def extract_tsd_centered_motifs_5p(
    tsv_path: str,
    genome_fasta: str,
    tsd_len: int = 5,
    flank_left: int = 10,
    flank_right: int = 10,
    min_mapq: Optional[int] = None,
    orient_to_plus: bool = True,
) -> pd.DataFrame:
    """
    DONOR-SIDE = 5' junction reads.

    ins0 in TSV is the *reference 5' end* of the read:
      strand '+' => ins0 corresponds to LEFT edge of the TSD copy adjacent to that junction
      strand '-' => ins0 corresponds to RIGHT edge of the TSD copy adjacent to that junction

    We compute tsd_start0 per row accordingly, then extract:
      [flank_left] + [TSD (tsd_len)] + [flank_right]

    sequences are oriented wrt top strand insertions, so insertions should be on the left hand side of the PAM

    """
    df = pd.read_csv(tsv_path, sep="\t")

    required = {"ref", "ins0", "strand"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"TSV missing required columns: {sorted(missing)}")

    if min_mapq is not None:
        if "mapq" not in df.columns:
            raise ValueError("min_mapq was set but TSV has no 'mapq' column.")
        df = df[df["mapq"] >= min_mapq].copy()

    # --- Compute true TSD start (0-based) for each read, strand-aware ---
    ins0 = df["ins0"].astype(int)
    strand = df["strand"].astype(str)

    # '+' : ins0 is left edge => tsd_start = ins0
    # '-' : ins0 is right edge => tsd_start = ins0 - (tsd_len - 1)
    df["tsd_start0"] = ins0.where(strand == "+", ins0 - (tsd_len - 1))

    # Drop impossible negative coordinates
    df = df[df["tsd_start0"] >= 0].copy()

    genome = hp.load_fasta_as_dict(genome_fasta)

    motifs = []
    keep = []

    for idx, row in df.iterrows():
        ref = str(row["ref"])
        if ref not in genome: #helpful if there are multiple sequences associated with a reference. since we are using a e.coli genome there is only ever one sequence.
            continue

        tsd_start0 = int(row["tsd_start0"])
        tsd_end_excl = tsd_start0 + tsd_len

        start = tsd_start0 - flank_left
        end_excl = tsd_end_excl + flank_right

        seq = genome[ref]
        if start < 0 or end_excl > len(seq):
            continue

        motif = seq[start:end_excl]  # flank_left + tsd_len + flank_right

        if orient_to_plus and row["strand"] == "-": #here a bottom strand INSERTION
            motif = hp.reverse_complement(motif)

        motifs.append(motif)
        keep.append(idx)

    out = df.loc[keep].copy()
    out["motif_seq"] = motifs
    out["motif_len"] = flank_left + tsd_len + flank_right
    out["flank_left"] = flank_left
    out["tsd_len"] = tsd_len
    out["flank_right"] = flank_right
    return out


def top_n_insertions_frequencies(insertions, N):
    """
    Return the top N most frequent insertion positions as frequencies.

    Parameters
    ----------
    insertions : list of int
        Insertion genomic positions.
    N : int
        Number of top positions to return.

    Returns
    -------
    list of tuple
        [(position, frequency), ...] sorted by descending frequency.
        Frequencies sum to <= 1.0.
    """
    if not insertions or N <= 0:
        return []

    total = len(insertions)
    counts = Counter(insertions)

    # Sort by frequency (desc), then position (asc)
    top = sorted(
        counts.items(),
        key=lambda x: (-x[1], x[0])
    )[:N]

    # Convert counts -> frequencies
    top_freqs = [(pos, count / total) for pos, count in top]

    return top_freqs

def plot_insertion_profile(tsv_file,gRNA,genome_sequence,fig_filename,bins=100):
    df = pd.read_csv(tsv_file, sep="\t")
    df = df[df["mapq"] >= 30]
    positions = df["ins0"]

    gmap=hp.map_guide_to_genome(gRNA,genome_sequence)
    # ---- plot ----
    fig, ax = plt.subplots(figsize=(10, 3))

    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'font.weight': 'bold',          # global default
    'axes.labelweight': 'bold',
    'axes.titleweight': 'bold',
})


    font = {'family': 'Arial', 'size': 24}

    ax.hist(
        positions,
        bins=bins)
    
    ax.set_xlabel("Genomic position (bp, 0-based)", fontdict=font)
    ax.set_ylabel("Insertion count", fontdict=font)
    ax.set_yscale("log")

    
    # Draw vertical lines marking guide region
    ax.axvspan(
        gmap[0]['start'],
        gmap[0]['end'],
        alpha=0.3,
        label=f"guideRNA ({gmap[0]['strand']})"
    )

    ax.annotate(
        "",                              # no text
        xy=(gmap[0]['start'], 0),             # arrow tip (on x-axis)
        xytext=(gmap[0]['start'], -0.05),      # arrow tail (slightly below axis)
        arrowprops=dict(
            arrowstyle="->",
            color="red",
            linewidth=4
        ),
        annotation_clip=False
    )
    window_bp = 100
    pct_ontarget = percent_within_distance(positions, gmap[0]['start'], window_bp)
    ndx = gmap[0]['start']
    #ax.set_title(f"Insertion site distribution: {pct_ontarget:.1f}% at site: {ndx:.1f} with window: {window_bp:.1f}")
    
    fig.tight_layout()
    fig.savefig(fig_filename)
    plt.show()

def characterize_insertion_from_position(tsv,PAM, guideRNA, genome_sequence, spacer_len=20): #here tsv is already loaded
    pos = tsv['ins0'] #input tsv has to contain a single row corresponding to a single insertion event, otherwise this will fail.
    #returns all possible insertions based on matches to the spacers within the expected window
    #further filtering needs to occur after callign this function to narrow it down to the 'right' off-target insertion. 
    window = 120
    ts_window = 10
    strand = tsv['strand']
    if strand == '+':
        segment_to_search = genome_sequence[pos:pos+window]
        all_spacers = hp.extract_spacers( segment_to_search, 'NGG', spacer_len, '5prime', search_strands='-')
        #for any given system the spacer_side and PAM will be the same regardless of +/-
    else:
        segment_to_search = genome_sequence[pos-window:pos]
        all_spacers = hp.extract_spacers( segment_to_search, 'NGG', spacer_len, '5prime', search_strands='+')
    

    alns = []
    pam_dna = []
    pam_start = []
    pam_end = []
    protospacer_dna = []
    proto_start = []
    proto_end = []
    int_dna = []
    intervene_start = []
    intervene_end = []
    ts_dna = []
    ts_start = []
    ts_end = []
    insertion_length = []
    aln_score = []
    grna_strand = []
    
    df = pd.DataFrame(columns=["PAM","protospacer","gRNA_strand","intervening","insertion_length","target_site",
                               "PAM_start","PAM_end","proto_start","proto_end","aln_score",
                               "intervene_start","intervene_end","ts_start","ts_end"])
    
    for ii in range(0,len(all_spacers)):
        pam_dna.append(all_spacers[ii].pam_dna)
        protospacer_dna.append(all_spacers[ii].spacer_top) #sequences returned are always top strand
        
        if strand == '+': #insertion on top strand
            aa=hp.local_align( guideRNA, all_spacers[ii].spacer_top )
            alns.append( aa )
            aln_score.append( aa.score )
            pam_start.append(all_spacers[ii].pam_start+pos)
            pam_end.append(all_spacers[ii].pam_end+pos)
            proto_start.append(all_spacers[ii].spacer_start+pos)
            proto_end.append(all_spacers[ii].spacer_end+pos)
            intstart = pos
            intend = pos+all_spacers[ii].pam_start
            grna_strand.append('-')
        else:
            aa=hp.local_align( guideRNA, hp.reverse_complement( all_spacers[ii].spacer_top ))
            alns.append( aa )
            aln_score.append( aa.score )
            pam_start.append(all_spacers[ii].pam_start+pos-window)
            pam_end.append(all_spacers[ii].pam_end+pos-window)
            proto_start.append(all_spacers[ii].spacer_start+pos-window)
            proto_end.append(all_spacers[ii].spacer_end+pos-window)
            intstart = pos - window + all_spacers[ii].pam_end
            intend = pos
            grna_strand.append('+')
        insertion_length.append(intend-intstart+1)
        intervene_start.append(intstart) 
        intervene_end.append(intend)
        int_dna.append(genome_sequence[intstart:intend])
        ts_start.append(pos-ts_window)
        ts_end.append(pos+ts_window)
        ts_dna.append(genome_sequence[pos-ts_window:pos+ts_window])

    df['PAM'] = pam_dna
    df['protospacer']=protospacer_dna
    df['intervening']=int_dna
    df['target_site']=ts_dna
    df['PAM_start']=pam_start
    df['PAM_end']=pam_end
    df['proto_start']=proto_start
    df['proto_end']=proto_end
    df['intervene_start']=intervene_start
    df['intervene_end']=intervene_end
    df['ts_start']=ts_start
    df['ts_end']=ts_end
    df['insertion_length']=insertion_length
    df['aln_score']=aln_score
    df['gRNA_strand']=grna_strand
        
    return df, alns

def characterize_insertion_from_spacer(tsv,PAM, guideRNA, genome_sequence):
    #given a spacer that maps uniquely to the genome, give the insertion features based on our previous profiling data and what we actually observe. So.. returned df corresponds to real data, but only will choose one example to compute sequence characteristics.
    #i.e. distance should be 71 bp +/- 5 on either side. Window is 10 so that target-site sequence features can be extracted
    #this only really works for native spacers. a new function needs to be written if spacer sequences cannot be uniquely mapped to the genome. 

    #from spacer find insertions that occur at the 'correct' distance: 
    distance = 71
    window = 10
    
    grna = hp.map_guide_to_genome(guideRNA,genome_sequence)
    if grna[0]['strand'] == '+':
        predicted_pos = grna[0]['end']+3+distance
    else:
        predicted_pos = grna[0]['start']-3-distance
    
    df = pd.DataFrame(columns=["PAM","protospacer","gRNA_strand","intervening","insertion_length","target_site",
                               "PAM_start","PAM_end","proto_start","proto_end","aln_score",
                               "intervene_start","intervene_end","ts_start","ts_end"])

    ontarget_insertions = between(tsv,predicted_pos-window,predicted_pos+window)
    
    if len(ontarget_insertions['ins0']) == 0:
        return df
    best_ontarget = statistics.mode(ontarget_insertions['ins0'])
    

    df['proto_start']= [ grna[0]['start'] ]
    df['proto_end']= [ grna[0]['end'] ]
    df['protospacer']= [ grna[0]['sequence'] ]
    if grna[0]['strand'] == '+':
        df['PAM_start'] = [ grna[0]['end'] ]
        df['PAM_end']= [ grna[0]['end']+3 ]
        df['PAM']= [ genome_sequence[ grna[0]['end'] : grna[0]['end']+3 ] ]
        df['intervene_start']= [ grna[0]['end']+3 ]
        df['intervene_end']= [ best_ontarget ]
        df['intervening']= [ genome_sequence[ grna[0]['end']+3 :  best_ontarget  ] ]
        df['target_site']= [ genome_sequence[best_ontarget-window:best_ontarget+window] ]
        df['insertion_length']= [ len(genome_sequence[ grna[0]['end']+3 :  best_ontarget  ]) ] #insertion length computed from end of PAM to insertion point

    else:
        df['PAM_start'] = [ grna[0]['start']-3 ]
        df['PAM_end']= [ grna[0]['start'] ]
        df['PAM']= [ genome_sequence[grna[0]['start']-3:grna[0]['start']] ] #always list top strand of sequences
        df['intervene_start']= [ best_ontarget ]
        df['intervene_end']= [ grna[0]['start']-3 ]
        df['intervening']= [ genome_sequence[best_ontarget:grna[0]['start']-3] ]
        df['target_site']= [ genome_sequence[best_ontarget-window:best_ontarget+window] ]
        df['insertion_length']=[ len( genome_sequence[best_ontarget:grna[0]['start']-3] ) ]

    df['ts_start']= [ best_ontarget-window ]
    df['ts_end']= [ best_ontarget+window ]
    df['gRNA_strand']= [ grna[0]['strand'] ]
    return df

