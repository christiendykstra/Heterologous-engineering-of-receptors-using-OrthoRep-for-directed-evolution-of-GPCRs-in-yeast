# HERO analysis

Analysis code for:

> **Heterologous engineering of receptors using OrthoRep (HERO) for directed evolution of GPCRs in yeast**
> Dykstra C.B., Bean B.D.M., Rousseau O., Araujo F., Masoud D., Liu C.C., Whiteway M., Martin V.J.J.
> *bioRxiv*, 2026. doi: TODO preprint DOI

HERO is a platform for evolving GPCRs using OrthoRep continuous directed
evolution in yeast. HERO strains are passaged in the presence of agonist, so
that variants which signal more effectively grow better and are enriched.
Mutants are then isolated and screened to identify improved signallers, which
can be developed into biosensors.

This repository contains the sequencing pipeline used to genotype evolved
populations, the per-run configuration for every sequencing run in the study,
and the analysis scripts that generated the results reported in the paper.

Raw nanopore reads are deposited in the SRA under accession
**PRJNA1533496**.

---

## Scope

This repository contains the sequencing pipeline and the analysis code that
generated results reported in the paper. Figure panels assembled directly in
GraphPad Prism from the supplementary data tables are not represented here; the
figure map below states which panels came from which source.

---

## Layout

```
.
├── pipeline/                     Nanopore demux → consensus → mutation calling
│   ├── pipeline.sh               entry point; takes a pooled FASTQ
│   ├── call_mutations.py         step 5: consensus → per-well substitution calls
│   ├── barcode_arrangement_fwd.toml
│   ├── barcode_arrangement_rev.toml
│   ├── metadata/                 barcodes, well map, reference sequence
│   ├── environment.yml           conda environment for the pipeline
│   └── example_HTS20/            final calls table from one complete run
│
├── analysis_scripts/
│   ├── dose_response_curves/     4PL fitting, EC50 / Hill extraction
│   ├── mutational_analysis/      frequency, burden, substitution matrices, Ts/Tv, gene map
│   ├── mutation_comparison/      two-condition statistics (Mann-Whitney, K-S, Poisson RR, dN/dS)
│   └── substitution_ranker/      activation-category dot plot
│
├── environment-analysis.yml      conda environment for the analysis scripts
├── CITATION.cff
└── LICENSE
```

Each folder under `analysis_scripts/` holds one script plus an `example_fig*/`
folder containing a real input and the outputs it produces, so every script can
be run and verified without the full dataset.

---

## Installation

```bash
git clone https://github.com/christiendykstra/Heterologous-engineering-of-receptors-using-OrthoRep-for-directed-evolution-of-GPCRs-in-yeast.git
cd Heterologous-engineering-of-receptors-using-OrthoRep-for-directed-evolution-of-GPCRs-in-yeast

# analysis scripts
conda env create -f environment-analysis.yml
conda activate heroine-analysis

# sequencing pipeline (Linux)
conda env create -f pipeline/environment.yml
conda activate heroine-pipeline
```

Dorado is not a conda package. Install the standalone binary from
<https://github.com/nanoporetech/dorado>; v0.9.6+49e25e9 was used here.

---

## Running the pipeline

```bash
cd pipeline
bash pipeline.sh /path/to/reads.fastq.gz [output_dir]
```

Reads are demultiplexed in two passes (forward barcode, then reverse barcode on
reverse-complemented reads), grouped into wells via `metadata/barcode_groups.csv`,
aligned to the reference with minimap2, polished with medaka, and called against
the reference CDS by `call_mutations.py`.

Output lands in `pipeline_out/`. The file downstream scripts consume is
`pipeline_out/mutations/mutations_clean.csv` — one row per well, with
`NT_substitutions`, `AA_substitutions_nonsynonymous`,
`AA_substitutions_synonymous` and their counts.

`pipeline/example_fig5e/mutations_clean.csv` is that table for run 8, the
worked example this folder is configured for.

Tuning: `THREADS`, `MIN_READS`, `MIN_LEN` and `MEDAKA_MODEL` can be set in the
environment. Medaka uses the GPU when torch is built with CUDA; CPU inference
gives identical consensus calls but is considerably slower.

---

## Pipeline configuration and deposited runs

`pipeline/` is a complete worked example: the scripts plus the exact
configuration used for run 8. Every file under `pipeline/metadata/` is specific
to that run.

`pipeline_runs/` holds the corresponding configuration for each sequencing run
deposited under SRA accession PRJNA1533496. To process a different run, replace
`pipeline/metadata/` with that run's `metadata/` folder, drop its FASTQ into
`pipeline/`, and run as above.

| Run | Receptor | Figures |
|---|---|---|
| `run1` | OPRM1 | Fig. S6 |
| `run2` | OPRM1 | Fig. 2c,d; Fig. S7 |
| `run3` (plate1, plate2) | OPRM1 | Fig. S8 |
| `run4` (plate1, plate2) | OPRM1 | Fig. 3c,d |
| `run5` (plate1, plate2) | OPRM1 | Fig. 4h,i |
| `run6` (plate1, plate2) | OPRD1 | Fig. S19 |
| `run7` | OPRD1 | Fig. S20, S21 |
| `run8` | OPRD1 | Fig. 5e; Fig. S20, S21 |

Runs split into `plate1` and `plate2` were sequenced separately because both
plates used the same barcode set, so each plate has its own FASTQ and its own
well map.

Each `metadata/` folder is self-contained:

| File | What it defines |
|---|---|
| `barcode_groups.csv` | barcode pair → well map for that plate |
| `fwd_barcodes.fasta` / `rev_barcodes.fasta` | barcode sequences in use |
| `barcode_arrangement_fwd.toml` / `_rev.toml` | barcode index ranges and the flanking primer masks |
| `reference.fasta` | alignment reference (minimap2, medaka) |
| `reference.csv` | reference CDS in `ref_seq_NT` (mutation calling) |

`call_mutations.py` derives its CDS search anchor from `ref_seq_NT`, so changing
the reference is sufficient to retarget it — no code edit needed.

---

## Running the analysis scripts

Every script takes its input as a command-line argument and writes outputs
beside it. Run any of them with `--help` for the full option list.

```bash
python analysis_scripts/dose_response_curves/doseresponse.py data.xlsx
python analysis_scripts/mutational_analysis/mutanal.py calls.xlsx
python analysis_scripts/mutation_comparison/mutation_comparison.py calls.xlsx
python analysis_scripts/substitution_ranker/substitution_ranker.py screen.xlsx
```

Input formats are documented in each script's docstring, and each
`example_fig*/` folder contains a working input.

### dose_response_curves
Fits four-parameter logistic dose-response curves; exports figures plus a CSV of
EC50, Hill slope, asymptotes, R² and confidence intervals. Input: one sheet,
column 0 = log10[ligand] in M, then N samples × 4 replicate columns.

### mutational_analysis
Mutation frequency distributions, burden histograms, position distributions, AA
and NT substitution matrices, a 4×4 NT matrix annotated with Ts/Tv, and a gene
map of nonsynonymous and synonymous density along the ORF.

### mutation_comparison
Compares per-clone mutation counts between two conditions: Mann-Whitney U with
rank-biserial effect size, Kolmogorov-Smirnov, exact-conditional Poisson rate
ratio with score-method 95% CI, Fisher exact on dN/dS, and the minimum rate
ratio detectable at 80% power.

### substitution_ranker
Places each carrier genotype of a substitution into the lowest agonist
concentration at which it clears the activation gate, then plots carriers per
category with marker size as carrier fraction and colour as mean µmax.

---

## Notes

- Figures are exported as SVG with editable text (`svg.fonttype = 'none'`) and
  were assembled into final panels in Illustrator. Differences between script
  output and published panels are layout only; no data were altered.
- `call_mutations.py` uses `Bio.pairwise2`, which is deprecated in Biopython.
  The pinned version in `pipeline/environment.yml` is required to reproduce
  these calls exactly.
- Amino-acid calls use reference-coordinate codon walking, so codons downstream
  of an indel stay in the reference reading frame. Indel-affected codons are
  skipped rather than mistranslated, and such wells are tagged
  `[indel_corrected]` in the full output table.
- Scripts set Arial as the default font where used. On systems without Arial,
  matplotlib falls back to a default sans-serif with a warning; figures render
  correctly otherwise.

---

## Citation

If you use this code, please cite the paper (above) and this repository:

```
https://doi.org/10.5281/zenodo.22903774
```

See `CITATION.cff` for machine-readable citation metadata.

## License

MIT — see `LICENSE`.

## Contact

Christien B. Dykstra — Centre for Applied Synthetic Biology, Concordia University
dykstrachristien@gmail.com

Questions and issues: please open a GitHub issue.
