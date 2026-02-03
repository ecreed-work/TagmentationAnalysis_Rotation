#!/usr/bin/env python3
'''
#     Generate insertion-centered sequence logos from junction mapping results
#     Usage example:

python reads_logo_from_manifest.py \
  --manifest-xlsx /Users/ecreed/Desktop/WorkingTagmentationAnalysis/test_output/filtered_tsv_manifest.xlsx \
  --outdir /Users/ecreed/Desktop/WorkingTagmentationAnalysis/test_output \
  --flank 100 \
  --logo-type bits \
  --logo-name W1DN1_analysis_logo+reads
  '''

import argparse
import os
from collections import Counter
import glob
import math

import pandas as pd
import matplotlib.pyplot as plt
import logomaker


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def load_fasta(fasta_path: str) -> dict:
    seqs = {}
    name = None
    chunks = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs[name] = "".join(chunks).upper()
    if not seqs:
        raise ValueError(f"No sequences found in FASTA: {fasta_path}")
    return seqs


def extract_window(ref_seq: str, ins0: int, flank: int) -> str:
    L = len(ref_seq)
    out = []
    for pos in range(ins0 - flank, ins0 + flank + 1):
        out.append(ref_seq[pos] if 0 <= pos < L else "N")
    return "".join(out)


def summarize_base_freq(windows: list[str]) -> pd.DataFrame:
    wlen = len(windows[0])
    counts = {b: [0] * wlen for b in ["A", "C", "G", "T", "N"]}

    for w in windows:
        for i, ch in enumerate(w):
            counts[ch if ch in counts else "N"][i] += 1

    df = pd.DataFrame(counts)
    df["total"] = df.sum(axis=1)

    for b in ["A", "C", "G", "T", "N"]:
        df[f"freq_{b}"] = df[b] / df["total"].replace(0, pd.NA)

    flank = wlen // 2
    df.insert(0, "pos_rel", range(-flank, flank + 1))
    return df


def bits_logo_matrix(freq_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in freq_df.iterrows():
        ps = [r[f"freq_{b}"] if pd.notna(r[f"freq_{b}"]) else 0.0 for b in ["A", "C", "G", "T"]]
        s = sum(ps)
        if s == 0:
            rows.append((r["pos_rel"], 0, 0, 0, 0))
            continue
        ps = [p / s for p in ps]
        H = -sum(p * math.log2(p) for p in ps if p > 0)
        IC = max(0.0, 2.0 - H)
        rows.append((r["pos_rel"], *(IC * p for p in ps)))

    return pd.DataFrame(rows, columns=["pos_rel", "A", "C", "G", "T"]).set_index("pos_rel")


def enriched_kmers(windows: list[str], k: int, max_kmers: int = 200) -> pd.DataFrame:
    c = Counter()
    for w in windows:
        for i in range(len(w) - k + 1):
            kmer = w[i:i + k]
            if "N" not in kmer:
                c[kmer] += 1
    return pd.DataFrame(c.most_common(max_kmers), columns=["kmer", "count"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest-xlsx", required=True)
    ap.add_argument("--sheet", default=0)
    ap.add_argument("--tsv-col", default="filtered_tsv")
    ap.add_argument("--ref-col", default="ref")
    ap.add_argument("--sample-col", default="samplename")
    ap.add_argument("--ins-col", default="ins0")
    ap.add_argument("--refname-col", default="ref")
    ap.add_argument("--strand-col", default="strand")
    ap.add_argument("--flank", type=int, default=50)
    ap.add_argument("--logo-type", choices=["frequency", "bits"], default="frequency")
    ap.add_argument("--logo-name", default="insertion_logo",
                    help="Base name used for plot title and output files")
    ap.add_argument("--kmer", type=int, default=6)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    manifest = pd.read_excel(args.manifest_xlsx, sheet_name=args.sheet)
    fasta_cache = {}

    windows = []
    reads_processed = 0

    for _, row in manifest.iterrows():
        tsv_path = str(row[args.tsv_col])
        ref_fasta = str(row[args.ref_col])

        if ref_fasta not in fasta_cache:
            fasta_cache[ref_fasta] = load_fasta(ref_fasta)
        ref_dict = fasta_cache[ref_fasta]

        for tsv in glob.glob(tsv_path):
            df = pd.read_csv(tsv, sep="\t")
            for _, r in df.iterrows():
                reads_processed += 1
                contig = r[args.refname_col]
                ins0 = int(r[args.ins_col])
                strand = r[args.strand_col]

                w = extract_window(ref_dict[contig], ins0, args.flank)
                if strand == "-":
                    w = revcomp(w)
                windows.append(w)

    if not windows:
        raise RuntimeError("No insertion windows extracted.")

    freq_df = summarize_base_freq(windows)
    freq_df.to_csv(os.path.join(args.outdir, "base_frequencies.tsv"),
                   sep="\t", index=False)

    if args.logo_type == "frequency":
        logo_df = freq_df[["pos_rel", "freq_A", "freq_C", "freq_G", "freq_T"]] \
            .rename(columns={"freq_A": "A", "freq_C": "C",
                             "freq_G": "G", "freq_T": "T"}) \
            .set_index("pos_rel")
        ylabel = "Base frequency"
    else:
        logo_df = bits_logo_matrix(freq_df)
        ylabel = "Information (bits)"

    png_out = os.path.join(args.outdir, f"{args.logo_name}_{args.logo_type}.png")
    pdf_out = os.path.join(args.outdir, f"{args.logo_name}_{args.logo_type}.pdf")

    plt.figure(figsize=(12, 3))
    logo = logomaker.Logo(logo_df)
    logo.ax.axvline(0, color="black", linestyle="--", linewidth=1)
    logo.ax.set_xlabel("Position relative to insertion site")
    logo.ax.set_ylabel(ylabel)
    logo.ax.set_title(
        f"{args.logo_name} — insertion-centered sequence logo ({args.logo_type})\n"
        f"Reads processed: {reads_processed} | Windows used: {len(windows)}"
    )
    plt.tight_layout()
    plt.savefig(png_out, dpi=300)
    plt.savefig(pdf_out)
    plt.close()

    enriched_kmers(windows, args.kmer).to_csv(
        os.path.join(args.outdir, f"enriched_{args.kmer}mers.tsv"),
        sep="\t", index=False
    )

    print(f"[ok] windows: {len(windows)}")
    print(f"[wrote] {png_out}")
    print(f"[wrote] {pdf_out}")


if __name__ == "__main__":
    main()