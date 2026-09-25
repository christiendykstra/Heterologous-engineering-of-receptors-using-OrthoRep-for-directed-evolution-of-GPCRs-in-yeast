#!/usr/bin/env python3
"""
call_mutations.py
=================

Copyright (c) 2026 Christien Dykstra
SPDX-License-Identifier: MIT

Call nucleotide and amino-acid mutations from per-well consensus sequences.

Step 5 of the nanopore amplicon pipeline. Takes the medaka consensus FASTA for
each well, locates the coding sequence, aligns it to the reference CDS, and
writes per-well substitution calls.

Amino-acid calls use reference-coordinate codon walking: alignment columns are
grouped into codons by REFERENCE position, so downstream codons stay in the
reference reading frame regardless of upstream indels. Codons containing a gap
in either sequence are skipped as indel-affected rather than mistranslated.


INPUT
-----
--outdir must contain the directories written by the pipeline:

    consensus/<well>_consensus.fasta    one consensus per well
    aligned/<well>.bam                  used only for the mapped-read count

--ref-csv is a one-row CSV describing the reference. Required column:

    ref_seq_NT      the reference coding sequence, starting at ATG

Other columns (ref_seq_name, ref_seq_alignment, ref_seq_coding, use_longest_ORF,
structure_file) are ignored by this script.

The consensus is searched for the first --anchor-len bases of ref_seq_NT to
locate the start of the CDS; both strands are tried. Wells where the anchor
isn't found are reported and skipped.


OUTPUT
------
Written to <outdir>/mutations/:

    mutations.csv           barcode, read_count, nt_mutations, aa_mutations,
                            n_nt_mutations, n_aa_mutations
                            NT calls include indels as ins<pos>_<base> /
                            del<pos>_<base>; AA calls carry an
                            [indel_corrected] tag when the well had any indel

    mutations_clean.csv     barcode, Count, NT_substitutions,
                            NT_substitutions_count,
                            AA_substitutions_nonsynonymous,
                            AA_substitutions_synonymous,
                            AA_substitutions_nonsynonymous_count
                            substitutions only, indels stripped. This is the
                            table the downstream analysis scripts consume.

Note on AA_substitutions_synonymous: entries are synonymous changes reported in
NUCLEOTIDE coordinates (e.g. C336T), since a synonymous change has no amino-acid
call to name. Downstream scripts use it as a count.


USAGE
-----
    python call_mutations.py --outdir pipeline_out --ref-csv metadata/reference.csv
    python call_mutations.py                     # uses ./pipeline_out and ./metadata/

Called automatically by pipeline.sh as step 5.

Requires: biopython
NOTE: this uses Bio.pairwise2, which is deprecated in Biopython and slated for
removal. Pin biopython <=1.88 (see environment.yml) to reproduce these calls;
switching to Bio.Align.PairwiseAligner may change alignments at indels.
"""

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

try:
    from Bio import SeqIO, pairwise2
    from Bio.Seq import Seq
except ImportError:
    sys.exit("ERROR: Biopython not found. Install with: "
             "conda install -c conda-forge biopython")

SCRIPT_DIR = Path(__file__).resolve().parent

# Alignment scoring: match, mismatch, gap-open, gap-extend
MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 2, -1, -2, -0.5

# How many bases from the start of the reference CDS are used to locate the
# coding sequence inside each consensus.
ANCHOR_LEN = 13


def parse_args():
    p = argparse.ArgumentParser(
        description="Call NT and AA mutations from per-well consensus sequences.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="See the module docstring for the expected directory layout.")
    p.add_argument("--outdir", default=os.environ.get("OUTDIR",
                   str(SCRIPT_DIR / "pipeline_out")),
                   help="Pipeline output dir containing consensus/ and aligned/")
    p.add_argument("--ref-csv", default=os.environ.get("REF_CSV",
                   str(SCRIPT_DIR / "metadata" / "reference.csv")),
                   help="One-row CSV with a ref_seq_NT column")
    p.add_argument("--anchor-len", type=int, default=ANCHOR_LEN,
                   help=f"Bases of the reference CDS used to locate the start "
                        f"codon in each consensus (default: {ANCHOR_LEN})")
    return p.parse_args()


def load_reference(ref_csv, anchor_len):
    """Read the reference CDS and derive the search anchor from it."""
    with open(ref_csv, newline="") as f:
        reader = csv.DictReader(f)
        row = next(reader, None)

    if row is None:
        sys.exit(f"ERROR: reference CSV is empty: {ref_csv}")
    if "ref_seq_NT" not in row:
        sys.exit(f"ERROR: reference CSV has no 'ref_seq_NT' column. "
                 f"Found: {list(row.keys())}")

    ref_coding = str(row["ref_seq_NT"]).strip().upper()
    if not ref_coding:
        sys.exit(f"ERROR: ref_seq_NT is empty in {ref_csv}")

    anchor = ref_coding[:anchor_len]
    ref_aa = str(Seq(ref_coding).translate(to_stop=True))
    print(f"  Reference CDS: {len(ref_coding)} nt → {len(ref_aa)} aa")
    print(f"  CDS anchor: {anchor}")
    return ref_coding, anchor


def mapped_read_count(bam_file):
    """Mapped reads in a BAM, or 0 if the file or samtools is unavailable."""
    if not bam_file.exists():
        return 0
    try:
        result = subprocess.run(
            ["samtools", "view", "-c", "-F", "4", str(bam_file)],
            capture_output=True, text=True)
    except FileNotFoundError:
        return 0
    return int(result.stdout.strip()) if result.returncode == 0 else 0


def find_coding(seq, anchor, ref_len):
    """Locate the CDS in a consensus by its anchor, trimmed to whole codons."""
    pos = seq.find(anchor)
    if pos == -1:
        return None
    coding = seq[pos:]
    trim = min(len(coding), ref_len)
    trim = trim - trim % 3
    return coding[:trim]


def call_well(record, ref_coding, anchor, barcode, read_count):
    """Align one consensus to the reference and return its mutation calls."""
    query_full = str(record.seq).upper()
    query_rc = str(Seq(query_full).reverse_complement())

    query_coding = find_coding(query_full, anchor, len(ref_coding))
    orientation = "fwd"
    if query_coding is None:
        query_coding = find_coding(query_rc, anchor, len(ref_coding))
        orientation = "rc"
        if query_coding is None:
            print(f"  WARNING: CDS anchor not found in {barcode} — skipping")
            return None

    print(f"  {barcode} ({read_count} reads): using {orientation} strand")

    trim = len(query_coding)
    alignments = pairwise2.align.globalms(
        ref_coding[:trim], query_coding,
        MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND,
        one_alignment_only=True)

    if not alignments:
        print(f"  WARNING: Alignment failed for {barcode}")
        return None

    aln_ref, aln_query = alignments[0][0], alignments[0][1]

    # --- NT mutations (substitutions + indels, 1-based reference coords) ---
    nt_mutations = []
    ref_pos = 0
    for r, q in zip(aln_ref, aln_query):
        if r != "-":
            ref_pos += 1
        if r == "-" and q != "-":
            nt_mutations.append(f"ins{ref_pos}_{q}")
        elif q == "-" and r != "-":
            nt_mutations.append(f"del{ref_pos}_{r}")
        elif r != q:
            nt_mutations.append(f"{r}{ref_pos}{q}")

    has_indels = any(m.startswith(("ins", "del")) for m in nt_mutations)

    # --- group alignment columns into codons by REFERENCE position ---
    # Every 3 reference bases = 1 codon. Insertions don't consume a reference
    # position, so they are recorded as a gap marker in the current codon,
    # which flags it as indel-affected. Downstream codons therefore stay in
    # the reference reading frame regardless of upstream indels.
    codon_ref, codon_query, codons = [], [], []
    for r, q in zip(aln_ref, aln_query):
        if r == "-":
            codon_query.append("-")   # query insertion; no ref position consumed
            continue
        codon_ref.append(r)
        codon_query.append(q)         # a base, or "-" for a deletion
        if len(codon_ref) == 3:
            codons.append(("".join(codon_ref), "".join(codon_query)))
            codon_ref, codon_query = [], []

    # --- AA substitutions and synonymous NT changes ---
    aa_mutations, syn_muts = [], []
    for i, (rc, qc) in enumerate(codons):
        if "-" in rc or "-" in qc:
            continue                  # indel-affected codon
        ref_aa_i = str(Seq(rc).translate())
        query_aa_i = str(Seq(qc).translate())
        if ref_aa_i != query_aa_i:
            aa_mutations.append(f"{ref_aa_i}{i + 1}{query_aa_i}")
        elif rc != qc:
            # codon changed, amino acid did not: report in NT coordinates
            for j in range(3):
                if rc[j] != qc[j]:
                    syn_muts.append(f"{rc[j]}{i * 3 + j + 1}{qc[j]}")

    nt_str = ", ".join(nt_mutations) if nt_mutations else "WT"
    aa_str = ", ".join(aa_mutations) if aa_mutations else "WT"
    if has_indels:
        aa_str += " [indel_corrected]"

    nt_subs_only = [m for m in nt_mutations if not m.startswith(("ins", "del"))]

    full = {
        "barcode": barcode,
        "read_count": read_count,
        "nt_mutations": nt_str,
        "aa_mutations": aa_str,
        "n_nt_mutations": len(nt_mutations),
        "n_aa_mutations": len(aa_mutations),
    }
    clean = {
        "barcode": barcode,
        "Count": read_count,
        "NT_substitutions": ", ".join(nt_subs_only) if nt_subs_only else "WT",
        "NT_substitutions_count": len(nt_subs_only),
        "AA_substitutions_nonsynonymous": ", ".join(aa_mutations) if aa_mutations else "WT",
        "AA_substitutions_synonymous": ", ".join(syn_muts),
        "AA_substitutions_nonsynonymous_count": len(aa_mutations),
    }

    print(f"    NT={nt_str[:80]}{'...' if len(nt_str) > 80 else ''}")
    print(f"    AA={aa_str[:80]}{'...' if len(aa_str) > 80 else ''}")
    return full, clean


def main():
    args = parse_args()

    outdir = Path(args.outdir).expanduser().resolve()
    ref_csv = Path(args.ref_csv).expanduser().resolve()

    if not ref_csv.is_file():
        sys.exit(f"ERROR: reference CSV not found: {ref_csv}")

    consensus_dir = outdir / "consensus"
    aligned_dir = outdir / "aligned"
    mutations_dir = outdir / "mutations"

    if not consensus_dir.is_dir():
        sys.exit(f"ERROR: no consensus directory in {outdir}. "
                 f"Run the pipeline through step 4 first.")
    mutations_dir.mkdir(parents=True, exist_ok=True)

    ref_coding, anchor = load_reference(ref_csv, args.anchor_len)

    fasta_files = sorted(consensus_dir.glob("*_consensus.fasta"))
    if not fasta_files:
        sys.exit(f"ERROR: no *_consensus.fasta files in {consensus_dir}")

    results, results_clean = [], []
    for fasta_file in fasta_files:
        barcode = fasta_file.stem.replace("_consensus", "")
        read_count = mapped_read_count(aligned_dir / f"{barcode}.bam")

        record = next(SeqIO.parse(fasta_file, "fasta"), None)
        if record is None:
            print(f"  WARNING: Empty consensus for {barcode}")
            continue

        called = call_well(record, ref_coding, anchor, barcode, read_count)
        if called is None:
            continue
        full, clean = called
        results.append(full)
        results_clean.append(clean)

    if not results:
        sys.exit("ERROR: no wells produced mutation calls.")

    fieldnames = ["barcode", "read_count", "nt_mutations", "aa_mutations",
                  "n_nt_mutations", "n_aa_mutations"]
    out_csv = mutations_dir / "mutations.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    clean_fieldnames = ["barcode", "Count", "NT_substitutions", "NT_substitutions_count",
                        "AA_substitutions_nonsynonymous", "AA_substitutions_synonymous",
                        "AA_substitutions_nonsynonymous_count"]
    out_csv_clean = mutations_dir / "mutations_clean.csv"
    with open(out_csv_clean, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=clean_fieldnames)
        writer.writeheader()
        writer.writerows(results_clean)

    print(f"\n  {len(results)} barcodes written to: {out_csv}")
    print(f"  {len(results_clean)} barcodes written to: {out_csv_clean} (indels stripped)")


if __name__ == "__main__":
    main()
