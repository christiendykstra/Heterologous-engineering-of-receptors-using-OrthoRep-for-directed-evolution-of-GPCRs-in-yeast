#!/usr/bin/env bash
# =============================================================================
# Nanopore amplicon pipeline: dual-barcode demux -> align -> consensus -> mutations
#
# Copyright (c) 2026 Christien Dykstra
# SPDX-License-Identifier: MIT
#
# Two-pass demultiplexing (Dorado FWD then REV) -> minimap2 -> medaka ->
# per-well mutation calling.
#
# USAGE
#   ./pipeline.sh <reads.fastq.gz> [output_dir]
#
#   ./pipeline.sh HTS20_pool.fastq.gz
#   ./pipeline.sh reads.fastq.gz /scratch/run1_out
#
# Raw reads are NOT included in this repository; download them from the SRA
# accession listed in the README. Everything else the pipeline needs lives in
# metadata/ and the two barcode arrangement TOMLs beside this script.
#
# REQUIRES
#   dorado, minimap2, samtools, medaka, python3 (+ biopython)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    echo "Usage: $(basename "$0") <reads.fastq.gz> [output_dir]" >&2
    exit 1
}

[[ $# -ge 1 ]] || usage
FASTQ="$1"
OUTDIR="${2:-$SCRIPT_DIR/pipeline_out}"

[[ -f "$FASTQ" ]] || { echo "ERROR: reads not found: $FASTQ" >&2; exit 1; }
FASTQ="$(cd "$(dirname "$FASTQ")" && pwd)/$(basename "$FASTQ")"

# ---------------------------------------------------------------------------
# CONFIG — metadata paths are relative to this script; tune the rest per run
# ---------------------------------------------------------------------------
FWD_BARCODE_FASTA="$SCRIPT_DIR/metadata/fwd_barcodes.fasta"
REV_BARCODE_FASTA="$SCRIPT_DIR/metadata/rev_barcodes.fasta"
FWD_TOML="$SCRIPT_DIR/metadata/barcode_arrangement_fwd.toml"
REV_TOML="$SCRIPT_DIR/metadata/barcode_arrangement_rev.toml"
BARCODE_GROUPS="$SCRIPT_DIR/metadata/barcode_groups.csv"
REFERENCE="$SCRIPT_DIR/metadata/reference.fasta"
REF_CSV="$SCRIPT_DIR/metadata/reference.csv"

THREADS="${THREADS:-8}"
MIN_READS="${MIN_READS:-10}"        # minimum reads to keep a well
MIN_LEN="${MIN_LEN:-1200}"          # minimum read length passed to REV demux
MEDAKA_MODEL="${MEDAKA_MODEL:-r1041_e82_400bps_sup_v5.2.0}"

for f in "$FWD_BARCODE_FASTA" "$REV_BARCODE_FASTA" "$FWD_TOML" "$REV_TOML" \
         "$BARCODE_GROUPS" "$REFERENCE" "$REF_CSV"; do
    [[ -f "$f" ]] || { echo "ERROR: missing required file: $f" >&2; exit 1; }
done

for tool in dorado minimap2 samtools medaka_consensus python3; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: '$tool' not found on PATH" >&2; exit 1; }
done

# ---------------------------------------------------------------------------
# SETUP — clean previous run, then create fresh directories
# ---------------------------------------------------------------------------
if [[ -z "$OUTDIR" || "$OUTDIR" == "/" ]]; then
    echo "ERROR: refusing to operate on output dir '$OUTDIR'" >&2; exit 1
fi
mkdir -p "$OUTDIR"
rm -rf "${OUTDIR:?}"/{demuxed_fwd,demuxed_paired,demuxed_wells,aligned,consensus,mutations}
mkdir -p "$OUTDIR"/{demuxed_fwd,demuxed_paired,demuxed_wells,aligned,consensus,mutations}

echo "============================================"
echo " Pipeline starting: $(date)"
echo " Reads:  $FASTQ"
echo " Output: $OUTDIR"
echo "============================================"

# ---------------------------------------------------------------------------
# STEP 1a: Demux by FWD barcode
# ---------------------------------------------------------------------------
echo ""
echo "[STEP 1a] Demultiplexing by FWD barcode..."

dorado demux \
    --output-dir "$OUTDIR/demuxed_fwd" \
    --kit-name fwd_barcode \
    --barcode-arrangement "$FWD_TOML" \
    --barcode-sequences "$FWD_BARCODE_FASTA" \
    --emit-fastq \
    --no-trim \
    "$FASTQ"

echo "  FWD demux results:"
for f in "$OUTDIR/demuxed_fwd/"*.fastq; do
    COUNT=$(grep -c "^@" "$f" || true)
    echo "    $(basename "$f"): $COUNT reads"
done

# ---------------------------------------------------------------------------
# STEP 1b: For each FWD barcode, demux by REV barcode
# ---------------------------------------------------------------------------
echo ""
echo "[STEP 1b] Demultiplexing by REV barcode..."

for FWD_FASTQ in "$OUTDIR/demuxed_fwd/"*.fastq; do
    FWD_NAME=$(basename "$FWD_FASTQ" .fastq)

    if [[ "$FWD_NAME" == *"unclassified"* ]]; then
        continue
    fi

    READ_COUNT=$(grep -c "^@" "$FWD_FASTQ" || true)
    if [ "$READ_COUNT" -lt 1 ]; then
        continue
    fi

    # Extract FWD barcode index
    FWD_IDX=$(echo "$FWD_NAME" | grep -oP 'barcode\K\d+')

    # Filter reads shorter than MIN_LEN before REV demux
    FILTERED_FASTQ="$OUTDIR/demuxed_fwd/${FWD_NAME}_filtered.fastq"
    awk -v min="$MIN_LEN" 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0} NR%4==0{if(length(s)>=min) print h"\n"s"\n"p"\n"$0}' \
        "$FWD_FASTQ" > "$FILTERED_FASTQ"
    FILT_COUNT=$(grep -c "^@" "$FILTERED_FASTQ" || true)
    echo "  FWD${FWD_IDX}: $READ_COUNT reads → $FILT_COUNT after length filter (≥${MIN_LEN}bp)"
    FWD_FASTQ="$FILTERED_FASTQ"
    READ_COUNT=$FILT_COUNT

    # Reverse-complement reads so REV barcode is at the 5' end for dorado
    RC_FASTQ="$OUTDIR/demuxed_fwd/${FWD_NAME}_rc.fastq"
    awk 'NR%4==1{print} NR%4==2{
        seq=$0; rc="";
        for(i=length(seq);i>=1;i--){
            c=substr(seq,i,1);
            if(c=="A")rc=rc"T"; else if(c=="T")rc=rc"A";
            else if(c=="G")rc=rc"C"; else if(c=="C")rc=rc"G";
            else rc=rc"N"
        } print rc
    } NR%4==3{print} NR%4==0{
        q=$0; rq="";
        for(i=length(q);i>=1;i--) rq=rq substr(q,i,1);
        print rq
    }' "$FWD_FASTQ" > "$RC_FASTQ"
    echo "  Reverse-complemented reads for REV demux"

    echo "  REV demuxing FWD barcode $FWD_IDX ($READ_COUNT reads)..."

    REV_OUTDIR="$OUTDIR/demuxed_paired/fwd${FWD_IDX}"
    mkdir -p "$REV_OUTDIR"

    dorado demux \
        --output-dir "$REV_OUTDIR" \
        --kit-name rev_barcode \
        --barcode-arrangement "$REV_TOML" \
        --barcode-sequences "$REV_BARCODE_FASTA" \
        --emit-fastq \
        --no-trim \
        "$RC_FASTQ" 2>/dev/null || true

    for f in "$REV_OUTDIR/"*.fastq; do
        COUNT=$(grep -c "^@" "$f" || true)
        REV_NAME=$(basename "$f" .fastq)
        if [[ "$REV_NAME" != *"unclassified"* ]] && [ "$COUNT" -gt 0 ]; then
            echo "    FWD${FWD_IDX} + $REV_NAME: $COUNT reads"
        fi
    done
done

# ---------------------------------------------------------------------------
# STEP 2: Resolve well IDs and combine paired FASTQ files
# ---------------------------------------------------------------------------
echo ""
echo "[STEP 2] Resolving well IDs from barcode_groups.csv..."

export OUTDIR BARCODE_GROUPS MIN_READS

# Temporarily relax strict mode for the Python heredoc
set +e
python3 << 'PYTHON_EOF'
import os
import csv
import re
from pathlib import Path

OUTDIR         = os.environ["OUTDIR"]
BARCODE_GROUPS = os.environ["BARCODE_GROUPS"]
MIN_READS      = int(os.environ["MIN_READS"])

# Load barcode groups: map (fwd_idx, rev_idx) -> well
bc_map = {}
with open(BARCODE_GROUPS, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        well = row["barcode_group"]
        fwd_idx = int(row["fwd"].replace("bc", ""))
        rev_idx = int(row["rvs"].replace("bc", ""))
        bc_map[(fwd_idx, rev_idx)] = well

paired_dir = Path(OUTDIR) / "demuxed_paired"
well_dir   = Path(OUTDIR) / "demuxed_wells"
well_dir.mkdir(exist_ok=True)

well_reads = {}

for fwd_dir in sorted(paired_dir.iterdir()):
    if not fwd_dir.is_dir():
        continue
    fwd_idx = int(fwd_dir.name.replace("fwd", ""))

    for fastq in sorted(fwd_dir.glob("*.fastq")):
        if "unclassified" in fastq.name:
            continue
        m = re.search(r'barcode(\d+)', fastq.name)
        if not m:
            continue
        rev_idx = int(m.group(1))

        key = (fwd_idx, rev_idx)
        if key not in bc_map:
            continue

        well = bc_map[key]
        count = sum(1 for line in open(fastq) if line.startswith("@"))

        if count < MIN_READS:
            print(f"  Skipping {well} (FWD{fwd_idx:02d}+REV{rev_idx:02d}): {count} reads < {MIN_READS}")
            continue

        well_reads.setdefault(well, []).append(fastq)
        print(f"  {well} (FWD{fwd_idx:02d}+REV{rev_idx:02d}): {count} reads")

# Concatenate reads for each well
for well, fastqs in sorted(well_reads.items()):
    out_fastq = well_dir / f"{well}.fastq"
    with open(out_fastq, "w") as out:
        for fq in fastqs:
            with open(fq) as inp:
                out.write(inp.read())
    total = sum(1 for line in open(out_fastq) if line.startswith("@"))
    print(f"  → Written {well}.fastq ({total} reads)")

PYTHON_EOF
set -e

# Check that step 2 produced output
WELL_COUNT=$(ls "$OUTDIR/demuxed_wells/"*.fastq 2>/dev/null | wc -l)
if [ "$WELL_COUNT" -eq 0 ]; then
    echo "  ERROR: No well FASTQ files produced by step 2 — check barcode_groups.csv"
    exit 1
fi
echo "  → $WELL_COUNT wells resolved"

# ---------------------------------------------------------------------------
# STEP 3: Align each well to reference with minimap2
# ---------------------------------------------------------------------------
echo ""
echo "[STEP 3] Aligning with minimap2..."

for FASTQ_WELL in "$OUTDIR/demuxed_wells/"*.fastq; do
    WELL=$(basename "$FASTQ_WELL" .fastq)
    READ_COUNT=$(grep -c "^@" "$FASTQ_WELL" || true)

    echo "  Aligning $WELL ($READ_COUNT reads)..."

    minimap2 \
        -ax map-ont \
        -t "$THREADS" \
        "$REFERENCE" \
        "$FASTQ_WELL" \
    | samtools sort -@ "$THREADS" -o "$OUTDIR/aligned/${WELL}.bam"

    samtools index "$OUTDIR/aligned/${WELL}.bam"

    MAPPED=$(samtools view -c -F 4 "$OUTDIR/aligned/${WELL}.bam")
    echo "    → $MAPPED reads mapped"
done

# ---------------------------------------------------------------------------
# STEP 4: Build consensus with medaka
# ---------------------------------------------------------------------------
echo ""
echo "[STEP 4] Building consensus with medaka..."

for BAM in "$OUTDIR/aligned/"*.bam; do
    WELL=$(basename "$BAM" .bam)

    MAPPED=$(samtools view -c -F 4 "$BAM")
    if [ "$MAPPED" -lt "$MIN_READS" ]; then
        echo "  Skipping $WELL — only $MAPPED mapped reads"
        continue
    fi

    echo "  Building consensus for $WELL..."

    medaka_consensus \
        -i "$BAM" \
        -d "$REFERENCE" \
        -o "$OUTDIR/consensus/$WELL" \
        -m "$MEDAKA_MODEL" \
        -t "$THREADS" \
        2>/dev/null

    if [ -f "$OUTDIR/consensus/$WELL/consensus.fasta" ]; then
        sed "s/^>.*/>$WELL/" \
            "$OUTDIR/consensus/$WELL/consensus.fasta" \
            > "$OUTDIR/consensus/${WELL}_consensus.fasta"
        echo "    → Consensus written"
    else
        echo "    WARNING: No consensus for $WELL"
    fi
done

# ---------------------------------------------------------------------------
# STEP 5: Call mutations
# ---------------------------------------------------------------------------
echo ""
echo "[STEP 5] Calling mutations..."

python3 "$SCRIPT_DIR/call_mutations.py" --outdir "$OUTDIR" --ref-csv "$REF_CSV"

echo ""
echo "============================================"
echo " Pipeline complete: $(date)"
echo " Output: $OUTDIR/mutations/"
echo "============================================"
