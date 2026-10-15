# Magibu App

Desktop application for Magibu: graphical interface and engine in a single program, distributed as a
Windows executable that needs no installation. **Research use only.**

## Download

Windows: download `Magibu-<version>-win.zip` from the [Releases](https://github.com/GianlucaMattei/MAGIBU/releases)
page, unzip it anywhere and run `Magibu.exe`. No installation and no administrator rights are needed: R and all
packages are included.

**Guides with screenshots: [express tutorial](docs/EXPRESS_TUTORIAL.md) (the essentials in 6 steps) and [full tutorial](docs/TUTORIAL.md).**

## What it does

- **Prediction**: projects new samples onto a reference model (default: `capper_v2`) and reports the closest classes.
  Input: a single file, all the files in a folder, or selected files. Output, in the chosen folder:
  - `<name>.distances.xlsx`: sheets `Summary` (3 nearest neighbours, closest groups, majority k-NN class),
    `Distances` (every class for every sample), `MeanDist_wide`
  - `<name>.results.plots.pdf`: local UMAP and bar plot for each sample
  - `Projection/<name>_metadata_combined.tsv`, `Tsv/` (converted samples), `<name>.magibu_result.rds`
    (to reopen the analysis in the app)
- **Global map**: interactive UMAP of the reference and of the projected samples.
- **Class zoom**: chosen classes and samples; the window grows to contain all the references of those classes and
  the samples. Shape = class, colour = subclass. Export to PDF/PNG.
- **Build model**: a folder with the training data (idat or Magibu tsv) and an annotation file (sample column,
  class column).

Interface languages: English, German, French, Spanish, Italian, Chinese (flags at the top right). Light and dark theme.

## Accepted input

| Input | Notes |
|---|---|
| Illumina idat (`*_Grn.idat` / `*_Red.idat`, also `.gz`) | EPIC v2, EPIC v1, 450k (array type detected or chosen) |
| Nanopore BAM | aligned to hg38, with MM/ML methylation tags (e.g. dorado 5mC or 5mC+5hmC), coordinate sorted, with the `.bai` index next to it |
| Entropy tsv (`*_sorted.entropy.file.tsv`) | output of `magibu-meth --entropy` or of modkit extract + ModkitResorter.sh (PoreMeth2); beta = column 5 |
| Magibu tsv | columns `Probe_ID`, `Beta`, `chr_hg38`, `pos_hg38` |

### Nanopore BAM

BAM files are read by `magibu-meth` (C, htslib; `app/bin/`, source in `build/native/prog/magibu_meth.c`), only in
the regions within 600 bp of the CpGs used by the chosen model, so a whole-genome BAM takes a few minutes.
No reference FASTA is needed: CpG sites are taken from the reads (CG in at least 20% of the reads).

Each read-level call is classified with the rules of the original Magibu nanopore workflow (ParseModkit.pl of
PoreMeth2), with p = (ML + 0.5) / 256:

| Option "Methylation" | Methylated | Unmethylated | Otherwise |
|---|---|---|---|
| 5mC only (default) | p5mC > 0.85 | p5mC + p5hmC < 0.3 | discarded |
| 5mC + 5hmC | p5mC + p5hmC > 0.85 | p5mC + p5hmC < 0.3 | discarded |

beta = methylated / (methylated + unmethylated), kept for CpGs with at least "Minimum coverage per CpG" calls
(default 1). With BAM files that carry only 5mC calls the two options give the same result.

The Summary reports `Model_CpG_coverage_pct`, the share of model CpGs measured in the sample (exact position or a
CpG within 500 bp). Below 50% the row is flagged: short runs (e.g. a few minutes of sequencing) do not cover enough
CpGs for a reliable result. BAM reading is available on Windows and Linux.

## magibu-meth on the command line

The same program can be used on its own (Linux binary in `app/bin/linux-x64/`, Windows in `app/bin/win-x64/`):

```bash
# whole genome entropy file for PoreMeth2 (replaces modkit extract full + ModkitResorter.sh)
magibu-meth --entropy --ref hg38.fasta -t 16 -o sample_sorted.entropy.file.tsv sample.bam
# beta only, in some regions, 5mC + 5hmC counted as methylation
magibu-meth --entropy --no-entropy --sum-5hmc --min-cov 1 -r regions.bed -t 8 -o sample.beta.tsv sample.bam
```

| Option | Meaning |
|---|---|
| `--ref FILE` | reference FASTA (with `.fai`); CpG sites from the reference instead of the reads |
| `-r FILE` | BED file: only these regions (default: chromosomes 1-22, X, Y) |
| `-t N` | threads (one chromosome per thread) |
| `--no-entropy` | beta only (entropy columns NA), no 3-CpG window |
| `--min-cov N` | minimum confident calls per CpG (default 3) |
| `--sum-5hmc` | methylated if p5mC + p5hmC > 0.85 |
| `--full-precision` | full precision numbers (default: 4 decimals) |

Output: tab separated, no header: chromosome, 0-based position of the C of the CpG (strands combined), entropy,
reads with the 3-CpG pattern, beta, coverage. On a 30x whole-genome BAM (134 GB) the entropy file was identical,
row by row, to the one of modkit extract + ModkitResorter.sh, in 6.5 minutes with 16 threads.

## Reproducibility

`uwot::umap_transform` is stochastic. By default the app behaves like the original pipeline (one batch, no seed):
projecting the same sample again moves it by about 0.2 UMAP units. `magibu_run_prediction(..., seed = 42)` projects
each sample alone with a fixed seed (same file, same result).

## Run from source (Linux)

Needs R (4.2 or later) with the packages in `build/packages.R`.

```bash
./run_local.sh 7654     # then open http://127.0.0.1:7654 (through an ssh tunnel if on a server)
```

## Folder structure

- `app/`: Shiny app. `app/R/magibu_lib.R` is the engine; `app/R/core/` holds copies of the Magibu scripts
  (two changes: sequential backend on Windows, RhpcBLASctl removed). `app/www/i18n.js`: translations.
- `models/capper_v2/`: reference model, reduced to what projection needs (20 MB instead of 1.6 GB, same projection).
- `electron/`: shell of the executable.
- `build/`: build scripts.

## Build the Windows executable

- On Linux (cross build, as for the released zip):
  ```bash
  build/win/prepare_win_r.sh       # downloads R 4.4.2 for Windows and the package binaries
  cd electron && npm install && npx electron-builder --win zip --x64 --config builder-win-cross.json
  ```
- On a Windows PC: `powershell -ExecutionPolicy Bypass -File build\build_windows.ps1`
- `magibu-meth`: `build/native/build_magibu_meth.sh` (gcc for Linux, zig cc for Windows; needs a static htslib 1.21
  built in `build/native/htslib-linux` and `build/native/htslib-win`). The compiled binaries are already in `app/bin/`.

## Licences

Third party software and licences: `build/THIRD_PARTY_LICENSES.txt` (R, Electron, R packages) and
`app/bin/LICENSES_magibu-meth.txt` (htslib, zlib, winpthreads). The `capper_v2` model is derived from the DNA
methylation reference cohort of Capper et al., Nature 2018 (GEO GSE90496). Research use only.
