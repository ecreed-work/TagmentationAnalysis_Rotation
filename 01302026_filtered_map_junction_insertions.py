#!/usr/bin/env python3
"""
execute with:
python ./01302026_filtered_map_junction_insertions.py \
  --fastq /Users/ecreed/Desktop/KelloggRotation/CodeProjects/tagmentation_inputdata/A14_S20_L001_R1_001.fastq \
  --ref /Users/ecreed/Desktop/KelloggRotation/CodeProjects/cJP003_assembly.fasta \
  --outdir /Users/ecreed/Desktop/KelloggRotation/CodeProjects/TagmentationAnalysis/01302026_output/test \
  --donor-seq TGGGTAAAGTCACA \
  --donor-side 5p \
  --min-read-count 100 \
  --min-mapq 30 \
  --qc-lowq-thresh 30 \ 
  --qc-max-lowq-frac 0.50 


Map donor–genome junction reads to a reference genome and call insertion coordinates per read.

Pipeline (unchanged mapping logic; QC only changes which reads are mapped):
  0) QC filter reads by Phred scores (FASTQ-level) -> writes qc-passed FASTQ
  1) Trim donor/transposon sequence using cutadapt
  2) Map to reference using bwa mem
  3) Sort/index BAM
  4) Call insertion site per read from BAM using CIGAR-aware coordinates
     + enforce perfect genome match (NM==0) with deviations allowed only as clipping (adaptor/junction)

Insertion coordinate definition (UNCHANGED):
  --donor-side 5p : insertion = reference 5' end of alignment
  --donor-side 3p : insertion = reference 3' end of alignment

Output TSV columns (same as original + appended QC metrics):
  read_id, ref, ins0, ins1, strand, mapq, cigar, aln_start0, aln_end0_excl,
  qc_read_len, qc_mean_q, qc_frac_q_lt20, qc_pct_q_ge30
-*
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
import gzip
import subprocess
from pathlib import Path
from typing import Dict, Optional, Tuple

import pysam


def run(cmd, *, shell=False):
    cmd_str = cmd if isinstance(cmd, str) else " ".join(map(str, cmd))
    print(f"[cmd] {cmd_str}")
    subprocess.run(cmd, check=True, shell=shell)


def open_maybe_gz(path: Path, mode: str):
    # FASTQ is typically text; use 'rt'/'wt'. Also supports .gz.
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def phred_scores(qual: str) -> list[int]:
    # Standard Illumina/Sanger Phred+33 encoding
    return [ord(c) - 33 for c in qual.rstrip("\n")]


def ensure_bwa_index(ref_fa: Path):
    # BWA index outputs: .amb .ann .bwt .pac .sa
    needed = [ref_fa.with_suffix(ref_fa.suffix + ext) for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]]
    if all(p.exists() for p in needed):
        return
    run(["bwa", "index", str(ref_fa)])


def qc_filter_fastq(
    fastq_in: Path,
    fastq_out: Path,
    *,
    lowq_thresh: int = 20,
    max_lowq_frac: float = 0.50,
) -> Dict[str, Tuple[int, float, float, float]]:
    """
    QC gate BEFORE mapping:
      drop if >= max_lowq_frac of bases have Q < lowq_thresh (default: >=50% bases Q<20)

    Writes QC-passing reads to fastq_out.

    Returns dict:
      read_id -> (read_len, mean_q, frac_q_lt20, pct_q_ge30)
    """
    metrics: Dict[str, Tuple[int, float, float, float]] = {}

    n_total = 0
    n_pass = 0

    with open_maybe_gz(fastq_in, "rt") as fin, open_maybe_gz(fastq_out, "wt") as fout:
        while True:
            h = fin.readline()
            if not h:
                break
            s = fin.readline()
            p = fin.readline()
            q = fin.readline()
            if not (s and p and q):
                raise ValueError("Truncated FASTQ record detected.")

            n_total += 1
            rid = h.strip().split()[0]
            if rid.startswith("@"):
                rid = rid[1:]

            qs = phred_scores(q)
            if not qs:
                continue

            read_len = len(qs)
            mean_q = sum(qs) / read_len
            n_lt20 = sum(1 for x in qs if x < lowq_thresh)
            frac_lt20 = n_lt20 / read_len
            n_ge30 = sum(1 for x in qs if x >= 30)
            pct_ge30 = 100.0 * (n_ge30 / read_len)

            metrics[rid] = (read_len, mean_q, frac_lt20, pct_ge30)

            if frac_lt20 >= max_lowq_frac:
                continue  # drop

            fout.write(h)
            fout.write(s)
            fout.write(p)
            fout.write(q)
            n_pass += 1

    print(f"[info] QC scanned reads: {n_total}")
    print(f"[info] QC passed reads: {n_pass}")
    if n_total:
        print(f"[info] QC pass rate: {n_pass / n_total:.4f}")
    print(f"[info] QC FASTQ: {fastq_out}")
    return metrics


def cutadapt_trim(
    fastq_in: Path,
    fastq_out: Path,
    donor_seq: str,
    donor_side: str,
    min_len: int,
    error_rate: float,
    cores: int,
):
    """Trim donor/transposon sequence from reads using cutadapt."""
    # cutadapt:
    # -g = 5' adapter, -a = 3' adapter
    if donor_side == "5p":
        adapter_arg = ["-g", donor_seq]
    else:
        adapter_arg = ["-a", donor_seq]

    cmd = [
        "cutadapt",
        *adapter_arg,
        "--discard-untrimmed",
        "-e", str(error_rate),
        "-m", str(min_len),
        "-j", str(cores),
        "-o", str(fastq_out),
        str(fastq_in),
    ]
    run(cmd)


def bwa_map(ref_fa: Path, fastq: Path, bam_out: Path, cores: int):
    """Map with bwa mem, then sort/index BAM."""
    sam_tmp = bam_out.with_suffix(".sam")

    # Write SAM via shell redirection (bwa writes SAM to stdout by default)
    run(f"bwa mem -t {cores} {ref_fa} {fastq} > {sam_tmp}", shell=True)

    run(["samtools", "sort", "-@", str(cores), "-o", str(bam_out), str(sam_tmp)])
    run(["samtools", "index", str(bam_out)])
    try:
        sam_tmp.unlink()
    except FileNotFoundError:
        pass


def ref_aln_span(read: pysam.AlignedSegment) -> Tuple[int, int]:
    return read.reference_start, read.reference_end


def insertion_coord(read: pysam.AlignedSegment, donor_side: str) -> Optional[int]:
    """Return insertion coordinate as 0-based integer, or None if cannot be called."""
    if read.is_unmapped:
        return None
    if read.reference_start is None or read.reference_end is None:
        return None
    if read.cigarstring is None or read.cigarstring == "*":
        return None

    start0, end0_excl = ref_aln_span(read)
    if end0_excl <= start0:
        return None

    forward = not read.is_reverse

    if donor_side == "5p":
        return start0 if forward else (end0_excl - 1)
    else:
        return (end0_excl - 1) if forward else start0


def is_perfect_genome_match_allow_clips(read: pysam.AlignedSegment) -> bool:
    """Require NM==0 and only {M, =, S, H} CIGAR ops."""
    if read.is_unmapped or read.cigartuples is None:
        return False

    try:
        nm = read.get_tag("NM")
    except KeyError:
        return False
    if nm != 0:
        return False

    allowed = {0, 4, 5, 7}  # M, S, H, =
    for op, _len in read.cigartuples:
        if op not in allowed:
            return False
    return True


def call_insertions(
    bam_path: Path,
    out_tsv: Path,
    *,
    donor_side: str,
    min_mapq: int,
    primary_only: bool,
    qc_metrics: Dict[str, Tuple[int, float, float, float]],
):
    bam = pysam.AlignmentFile(str(bam_path), "rb")

    n_total = 0
    n_kept = 0
    n_not_perfect = 0
    n_no_qc = 0

    with open(out_tsv, "w") as f:
        f.write("\t".join([
            "read_id", "ref", "ins0", "ins1", "strand", "mapq", "cigar", "aln_start0", "aln_end0_excl",
            "qc_read_len", "qc_mean_q", "qc_frac_q_lt20", "qc_pct_q_ge30",
        ]) + "\n")

        for read in bam.fetch(until_eof=True):
            n_total += 1

            if read.is_unmapped:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if primary_only and (read.is_secondary or read.is_supplementary):
                continue

            # Enforce perfect genomic match; allow only clipping for adaptor/junction
            if not is_perfect_genome_match_allow_clips(read):
                n_not_perfect += 1
                continue

            ins0 = insertion_coord(read, donor_side=donor_side)
            if ins0 is None:
                continue

            ref = bam.get_reference_name(read.reference_id)
            strand = "-" if read.is_reverse else "+"
            start0, end0_excl = ref_aln_span(read)

            rid = read.query_name
            qc = qc_metrics.get(rid)
            if qc is None:
                # Should be rare if BAM originated from QC-passed FASTQ; keep row but mark missing
                n_no_qc += 1
                qc_read_len, qc_mean_q, qc_frac_lt20, qc_pct_ge30 = ("NA", "NA", "NA", "NA")
            else:
                qc_read_len, qc_mean_q, qc_frac_lt20, qc_pct_ge30 = qc
                qc_mean_q = f"{qc_mean_q:.3f}"
                qc_frac_lt20 = f"{qc_frac_lt20:.4f}"
                qc_pct_ge30 = f"{qc_pct_ge30:.2f}"

            f.write("\t".join(map(str, [
                rid, ref, ins0, ins0 + 1, strand, read.mapping_quality,
                read.cigarstring, start0, end0_excl,
                qc_read_len, qc_mean_q, qc_frac_lt20, qc_pct_ge30,
            ])) + "\n")
            n_kept += 1

    bam.close()
    print(f"[info] BAM reads scanned: {n_total}")
    print(f"[info] Insertions written: {n_kept}")
    print(f"[info] Dropped (not perfect genome match): {n_not_perfect}")
    print(f"[info] Missing QC metrics in dict: {n_no_qc}")
    print(f"[info] Output TSV: {out_tsv}")


def filter_min_read_count_sites(
    tsv_in: Path,
    tsv_out: Path,
    *,
    min_read_count: int,
) -> tuple[int, int, int]:
    """Filter an insertion TSV by minimum reads per insertion *site*.

    A site is defined as (ref, ins0) and support is the number of TSV rows (reads) at that site.

    Returns: (n_rows_in, n_rows_out, n_sites_kept)
    """
    if min_read_count <= 0:
        raise ValueError("min_read_count must be > 0 to use this filter")

    site_counts: Counter[tuple[str, int]] = Counter()
    n_in = 0

    with open(tsv_in, "r", newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty TSV or missing header: {tsv_in}")
        if "ref" not in reader.fieldnames or "ins0" not in reader.fieldnames:
            raise ValueError(f"TSV missing required columns 'ref' and/or 'ins0': {tsv_in}")

        for row in reader:
            n_in += 1
            ref = row.get("ref")
            ins0_s = row.get("ins0")
            if not ref or ins0_s in (None, ""):
                continue
            try:
                ins0 = int(ins0_s)
            except ValueError:
                continue
            site_counts[(ref, ins0)] += 1

    kept_sites = {k for k, c in site_counts.items() if c >= min_read_count}

    n_out = 0
    with open(tsv_in, "r", newline="") as fin, open(tsv_out, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty TSV or missing header: {tsv_in}")
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for row in reader:
            ref = row.get("ref")
            ins0_s = row.get("ins0")
            if not ref or ins0_s in (None, ""):
                continue
            try:
                ins0 = int(ins0_s)
            except ValueError:
                continue
            if (ref, ins0) in kept_sites:
                writer.writerow(row)
                n_out += 1

    return n_in, n_out, len(kept_sites)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True, type=Path, help="Merged FASTQ (single-end).")
    ap.add_argument("--ref", required=True, type=Path, help="Reference genome FASTA (E. coli).")
    ap.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    ap.add_argument("--donor-seq", required=True, help="Donor/transposon sequence to trim (junction-proximal).")
    ap.add_argument("--donor-side", choices=["5p", "3p"], default="5p",
                    help="Where donor sequence is in the read before trimming: 5p or 3p.")

    ap.add_argument("--min-len", type=int, default=20, help="Minimum length of trimmed genomic read to keep.")
    ap.add_argument("--error-rate", type=float, default=0.1, help="cutadapt error rate (e.g. 0.1).")
    ap.add_argument("--cores", type=int, default=8, help="Threads for cutadapt/bwa/samtools.")

    ap.add_argument("--min-mapq", type=int, default=20, help="Minimum MAPQ to report insertion.")
    ap.add_argument("--primary-only", action="store_true", help="Drop secondary/supplementary alignments.")

    ap.add_argument("--qc-lowq-thresh", type=int, default=20, help="Low-Q threshold for QC filter (default 20).")
    ap.add_argument("--qc-max-lowq-frac", type=float, default=0.50,
                    help="QC drop if frac(Q<qc-lowq-thresh) >= this. Default 0.50.")

    # Analysis-layer filter: minimum reads supporting an insertion site (ref, ins0)
    ap.add_argument("--min-read-count", type=int, default=0,
                    help=("If >0, filter the output TSV to keep only insertion sites with at least this many reads. "
                          "If no sites pass, print \"remove from analysis\" and do not export a filtered TSV."))

    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    fastq = args.fastq
    ref = args.ref

    ensure_bwa_index(ref)

    # 0) QC filter
    qc_fastq = args.outdir / f"{fastq.stem}.qc.fastq"
    qc_metrics = qc_filter_fastq(
        fastq,
        qc_fastq,
        lowq_thresh=args.qc_lowq_thresh,
        max_lowq_frac=args.qc_max_lowq_frac,
    )

    # 1) Trim donor
    trimmed_fastq = args.outdir / f"{fastq.stem}.trimmed.fastq"
    cutadapt_trim(
        qc_fastq,
        trimmed_fastq,
        donor_seq=args.donor_seq,
        donor_side=args.donor_side,
        min_len=args.min_len,
        error_rate=args.error_rate,
        cores=args.cores,
    )

    # 2-3) Map + sort/index
    bam = args.outdir / f"{fastq.stem}.sorted.bam"
    bwa_map(ref, trimmed_fastq, bam, args.cores)

    # 4) Call insertion per read
    out_tsv = args.outdir / f"{fastq.stem}.insertions.{args.donor_side}.tsv"
    call_insertions(
        bam_path=bam,
        out_tsv=out_tsv,
        donor_side=args.donor_side,
        min_mapq=args.min_mapq,
        primary_only=args.primary_only,
        qc_metrics=qc_metrics,
    )

    # 5) Optional analysis-layer filter: minimum reads per insertion site (ref, ins0)
    if args.min_read_count and args.min_read_count > 0:
        filtered_tsv = args.outdir / f"filtered_{out_tsv.name}"
        n_in, n_out, n_sites = filter_min_read_count_sites(
            out_tsv, filtered_tsv, min_read_count=args.min_read_count
        )
        print(f"[info] Min-read-count filter (site-level): min_read_count={args.min_read_count}")
        print(f"[info] Input rows: {n_in}")
        print(f"[info] Sites kept: {n_sites}")
        print(f"[info] Output rows: {n_out}")

        if n_out == 0:
            print("remove from analysis")
            (args.outdir / f"{fastq.stem}.REMOVE_FROM_ANALYSIS.txt").write_text(
                f"remove from analysis\nmin_read_count={args.min_read_count}\n"
            )
            return

        print(f"[info] Filtered TSV: {filtered_tsv}")


if __name__ == "__main__":
    main()