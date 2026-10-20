# Changelog

All notable changes to MagCluster are documented here.

## [1.0.1] - 2026-10-20

Bug-fix release. MGC detection (gene prediction, HMM search, clustering rules and
thresholds) is unchanged. With default options, the result tables, protein FASTA and
GenBank slices are byte-identical to 1.0.0 for the same genome ID (see "Changed" for
how genome IDs are now derived).

### Fixed
- `-s/--score` had no effect: the per-profile `trusted_cutoff` (25.0) in `rules.yaml`
  always took precedence. An explicit `-s` value now applies to all profiles (including values
  below 25). Without `-s`, the per-profile cutoffs are used as before.
- Genome-name collisions: the genome ID was the file name up to the first `.`, so
  `A/genomic.gbff` and `B/genomic.gbff` (NCBI Datasets layout) were merged into one MGC.
  `strainX.v1.fna` and `strainX.v2.fna` also overwrote each other's genes. Genome IDs are now
  unique (see "Changed"). Hits and genes are kept per genome, so identical gene IDs in
  different genomes no longer mix.
- A directory containing a non-genome file (e.g. `README.txt`) crashed the whole batch and
  wrote no output. Only files with genome extensions are now read from directories; other
  files are skipped with a warning. A file that cannot be parsed is reported and skipped.
  Each genome's outputs are written as soon as it is processed.
- File names with spaces crashed GenBank/HTML writing (`Invalid whitespace ... LOCUS`).
  Genome IDs are now sanitised; the original path is recorded in the run log.
- Dependency floors: `typer>=0.16.0` and `pyhmmer>=0.10.1` (0.10.0 is yanked). typer
  <= 0.15.3 crashes (`TypeError ... make_metavar()`) together with click >= 8.2, which pip and
  conda install by default. typer 0.15.4 avoids this only through its own `click<8.2` pin.
  0.16.0 is the first release that works with every click version, and it was tested with
  click 8.0.0, 8.1.8, 8.2.0 and 8.5.0. The unused `pandas` and `pyfamsa` runtime
  dependencies were removed; `pyfamsa` is now the optional extra `magcluster[build-db]`,
  needed only by `scripts/build_database.py`.
- `--force` now removes the outputs of the previous run in that folder (tables, FASTA,
  `MGC_align.html`, `clusters_gbk/*.gbk`, `figures/*.svg`, `magcluster.log`) before
  writing. Other files are kept.
- A Clinker failure (e.g. a bad `--clinker-opts`) is now reported as an error with exit
  code 1. The run no longer claims that `MGC_align.html` was generated.
- Numeric options are validated with clear errors (exit code 2): `--threads >= 0`
  (values above the available cores are capped with a warning), `--evalue >= 0`,
  `--max-gap >= 0`, `--threshold >= 1`, `--flank >= 0`, `0 <= --identity <= 1`. Before,
  `-un -2` raised an IndexError and the other values were silently accepted.
- Clear one-line errors instead of Python tracebacks: empty, text or binary files given as
  FASTA; GFF3 input (not supported, and no longer documented as supported); unbalanced quotes
  in `--clinker-opts`; `-o` pointing to an existing file.
- `--single` on less than 20 kb of sequence crashed in Pyrodigal training. It now falls back
  to metagenome mode with a warning.
- MGCs on a complete circular replicon whose subclusters lie on both sides of the origin
  were reported with an almost genome-wide range (e.g. `5378-4365795`). They are now
  reported as the shortest arc (`start > end`, marked "spans origin"). Which genes belong to
  the MGC is unchanged.
- The same input listed twice (or a file plus its folder) is processed once. Before, its
  rows appeared twice in the table.
- `-h` now works as an alias of `--help`.
- `--clinker-help` showed "alignment scores ." because Rich swallowed `[default: 2]`.
- matplotlib and pyGenomeViz are imported only when `--static-plot` is used, and their font
  cache / font_manager messages are silenced. `--help`, `--version` and `--no-plot` runs
  print nothing from matplotlib. The root logger is no longer configured by MagCluster.
- Duplicate "(default: 3)" text removed from `--help`.
- `--meta/--single` now prints a notice when it is ignored because the input already has
  CDS annotations.
- Relative output paths in `--clinker-opts` (e.g. `-mo matrix.tsv`) are written inside the
  output directory instead of the current working directory.

### Added
- Run log: every `magcluster run` writes `magcluster.log` to the output directory
  (overwritten on re-runs). `--log FILE` writes it elsewhere. The log records the timestamp,
  MagCluster/Python/dependency versions, the full command line, all effective parameters,
  the genome ID of each input, per-genome progress, timings and results (including subcluster
  coordinates and families), and every warning or error (including skipped files). Console
  output is unchanged apart from a final "Run Log:" line.
- Python 3.13 in CI and package classifiers.
- Regression tests for all of the above, using a 90 kb MSR-1 magnetosome-island test genome
  (`tests/data/`).

### Changed
- Genome ID = file name minus `.gz` and one known genome extension, e.g.
  `GCF_000513295.1.fna` -> `GCF_000513295.1` (was `GCF_000513295`) and `strainX.v1.fna` ->
  `strainX.v1`. Whitespace and characters other than letters, digits, `.`, `-` and `_` become
  `_`. When two inputs give the same ID (compared case-insensitively), later ones get `_2`,
  `_3`, ... and a warning is printed. **Output file names, the `genome`/`cluster_id` columns and
  Pyrodigal gene IDs (`{genome}_00001`) change for inputs whose names contain extra dots or
  special characters.**
- Exit code 1 if any input could not be processed or Clinker failed. The results of the other
  genomes are still written.
- Directory inputs are processed in sorted file-name order.
- The `clinker` executable next to the running Python interpreter is preferred over `PATH`.

### Documentation
- README: install with `pip install magcluster` or
  `conda install -c conda-forge -c bioconda magcluster`. The logo now uses an absolute URL so
  it renders on PyPI. Option aliases (`--min-families`, `--unmag-num`, `-h`), all outputs
  (including `figures/` and `magcluster.log`) and the `-s` semantics are documented.
- README: corrected the MGC rule. An MGC needs at least `--threshold` (default 3) distinct
  magnetosome families, and **any non-Feo family counts**. It was described as "3 of the 10
  core families". The `clustering_params` and per-profile `evalue_cutoff` keys in `rules.yaml`
  are documented as informational (not read). They are left unchanged to avoid any behaviour
  change.

## [1.0.0] - 2026-10-19
- Prokka/BLAST-free rewrite: Pyrodigal gene prediction, PyHMMER search with 92 curated
  magnetosome HMM profiles, synteny-aware subcluster slicing, Clinker visualisation.
