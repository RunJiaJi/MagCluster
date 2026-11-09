# Changelog

All notable changes to MagCluster are documented here.

## [1.1.3] - 2026-11-09

HMM / sequence coverage reporting and `--min-hmm-cov` significance filtering release.

### Added
- **HMM and protein sequence coverage in `{genome}_mgc_genes.tsv`**:
  - Added `hmm_cov` (fraction of the HMM profile covered by aligned domain intervals, `0.00–1.00`)
    and `seq_cov` (fraction of the target protein covered by aligned domain intervals, `0.00–1.00`)
    to `MagnetosomeHit` and per-genome `{genome}_mgc_genes.tsv` tables, while keeping
    `{genome}_magnetosome_proteins.faa` headers concise.
- **Minimum HMM coverage threshold (`--min-hmm-cov`, default `0.5`)**:
  - Added `--min-hmm-cov` (`0.0–1.0`, default `0.5`) to `magcluster run` alongside `-e/--evalue`
    (`1e-5`) and `-s/--score` (`25.0`) to filter out spurious partial-motif matches (e.g. generic
    TPR/ATPase/protease short motif hits covering < 50% of the HMM model) and accurately assign
    single-domain proteins to domain profiles (such as `MamO_domain` vs `MamO_full`).

## [1.1.2] - 2026-11-09

Concise, non-redundant per-genome output formatting release. MGC calling rules, HMM profiles and
thresholds are unchanged.

### Changed
- **Concise, non-redundant per-genome TSV and FASTA outputs**:
  - Removed the redundant `genome` column from per-genome `{genome}_mgc_genes.tsv` tables
    (whose file names already carry the full genome name), shortened `cluster_id` inside
    per-genome tables and FASTA headers from `{genome}_MGC` to `MGC` (or `unclustered`), and
    dropped `[genome]` from `{genome}_magnetosome_proteins.faa` headers.
  - Pyrodigal-predicted gene IDs now use concise `gene_00001` numbering rather than repeating
    the full `{genome}_00001` prefix across `.tsv`, `.faa`, and `.gbk` `/locus_tag` qualifiers.

## [1.1.1] - 2026-11-07

Bug-fix and UX polish release following real-world user testing of v1.1.0. MGC calling rules,
HMM profiles and thresholds are unchanged.

### Fixed
- **Protein FASTA files with `.fasta`/`.fa` extensions**: files containing amino acid
  sequences (rather than DNA/RNA contigs) under `.fasta`, `.fa`, `.fna` or extensionless
  names are now automatically detected and parsed as protein FASTA with a notice, instead of
  being fed to Pyrodigal as nucleotide contigs and silently reporting `0` hits.
- **Duplicate sequence / locus IDs**: FASTA headers that share the same first whitespace
  token (e.g. `>genbankID| Magnetosome protein ... AGG16234.1 [FH-1]`) as well as duplicate
  contig IDs or duplicate GenBank `locus_tag` qualifiers are now automatically disambiguated
  and deduplicated with a warning, preventing hits from overwriting one another.
- **GenBank `Partial codon` BiopythonWarning**: `parse_genbank()` now skips untranslated
  `/pseudo` pseudogene CDS features, honours `codon_start` / `transl_table`, and trims
  trailing partial codons before *in silico* translation so NCBI `.gbff` files never trigger
  `BiopythonWarning: Partial codon, len(sequence) not a multiple of three`.
- **Multi-subcluster coordinate display & `_locN` numbering**:
  - `MGCCluster.format_coordinates()` now lists individual subcluster intervals whenever an
    MGC has multiple discrete subclusters—including multiple subclusters on the same
    chromosome (such as AMB-1's magnetosome island and magnetosome islet)—instead of
    reporting a chromosome-wide bounding span.
  - Fixed multi-contig `subcluster_id` numbering in `ClusterEngine.cluster_genome()` so
    subclusters across different contigs are numbered sequentially (`_loc1`, `_loc2`, ...) in
    `magcluster.log`, matching `.gbk` slices.
- **GenBank header consistency & Clinker locus labels**:
  - Genomes carrying a standard INSDC/RefSeq accession (e.g. `GCF_000009985.1_ASM998v1_genomic`)
    now consistently use that accession in `.gbk` `ACCESSION`/`VERSION`/`LOCUS` headers even
    when the full file name is under 68 characters; single-locus `DEFINITION` lines no longer
    repeat the word `cluster`.
  - `customize_clinker_html()` now populates untruncated, unique `ACCESSION` locus names
    (e.g. `GCF_000009985.1_loc1`..`loc4`, `MSR-1_genome_loc1`..`loc4`) in Clinker's HTML
    JSON instead of 16-character truncated `LOCUS` prefixes.
- **Rich markup stripping on batch progress lines**: escaped batch progress labels in
  `magcluster run` and `magcluster clinker` (`(batch 1/2, 50 clusters: first .. last)`) so
  Rich no longer strips `[batch ...]` as an unknown style tag.
- **No stray `magcluster_clinker.log` when `-o` targets an `.html` file**: passing a direct
  HTML file path (`magcluster clinker ... -o figure.html`) no longer drops an unrequested
  `magcluster_clinker.log` into the current working directory unless `--log` is explicitly
  provided (directory targets `-o <outdir>` continue writing `<outdir>/magcluster_clinker.log`).

### Added
- **Batch membership manifest (`MGC_align_batches.tsv`)**: whenever Clinker splits clusters
  into multiple HTML files (`n_batches > 1`), a companion TSV mapping each `cluster_id` to
  its `batch` number and `html_file` is written automatically (and cleaned up by `--force`).
- **Actionable Clinker diagnostics**:
  - Passing `-l/--list` with bare cluster/genome IDs without an input directory now hints to
    supply the GenBank/results directory alongside `-l`.
  - Passing an unsliced whole-genome GenBank file (`> 500` CDSs) to `magcluster clinker`
    emits a warning recommending `magcluster run` to slice MGC loci first.
  - When a generated Clinker HTML exceeds 30 MB due to dense pairwise alignments, a tip
    suggests lowering `--batch-size` (e.g. `-b 20`) or passing `--clinker-opts "-na"`.

## [1.1.0] - 2026-11-06

Feature and usability release for repository-scale MGC screening and standalone Clinker
visualization. MGC detection rules, HMM search parameters and output contents for normal
short genome names remain identical to 1.0.1.

### Added
- Standalone `magcluster clinker` subcommand: generate interactive Clinker HTML figures
  directly from sliced MGC GenBank files (`clusters_gbk/*.gbk`, `.gbff`, `.gbk.gz`),
  directories (including a `magcluster run` results folder, where `clusters_gbk/` is
  auto-detected), or a custom ordered list (`-l/--list FILE`). MagCluster's 8-family color
  palette and domain labels (`MamE_Cter`, `MamQ_greigite`, stripped `_full` suffixes) are
  restored automatically from embedded GenBank `/note` metadata without needing the original
  genomes or TSV tables; external third-party GenBank files are annotated on the fly with the
  built-in 92 magnetosome HMMs.
- Automatic Clinker batching (`-b/--batch-size`, default `50`; `0` disables splitting) in
  both `magcluster run` and `magcluster clinker`: when more than `--batch-size` MGC clusters
  are plotted, clusters are partitioned sequentially in input file order into
  `MGC_align_batch01.html`, `MGC_align_batch02.html`, ... so pairwise alignments stay fast
  and browser SVG rendering remains responsive. Runs with at most `--batch-size` MGC clusters
  continue to write a single `MGC_align.html`.
- Cohort summary table (`mgc_summary.tsv`) written to the output directory on every
  `magcluster run`, with one row per input genome reporting detection status, total genes,
  total and clustered HMM hits, cluster ID, subcluster count, coordinates, distinct family
  count, predicted mineral type (`Magnetite`, `Greigite`, `Putative/Mixed`), flags
  (`Contig_Edge`, `Spans_Origin`), family list, and sliced `.gbk` path.

### Fixed
- Long genome names (> 68 characters with `_MGC` or `_locN`) no longer trigger
  `BiopythonWarning: Annotation '...' too long for 'ACCESSION'/'VERSION' line` during GenBank
  export. While the `.gbk` file name (and thus Clinker's displayed cluster title) stays
  unchanged, internal GenBank `ACCESSION`/`VERSION`/`LOCUS` headers now retain only the
  extracted assembly/sequence accession (e.g. `GCA_002753255.1_loc1`, falling back to safe
  truncation when no standard accession is present), and `LOCUS` names preserve `_locN`
  suffixes across multi-locus slices. Full genome names remain recorded in `DEFINITION`
  (soft-wrapped across lines when exceeding 68 characters) and `/organism`.
- Multi-contig MGCs were displayed in console/log summaries with global `min(start)-max(end)`
  appended to a comma-joined contig list (e.g. `ctg1,ctg2,ctg3:3-7681`). Coordinates are now
  formatted per subcluster (`ctg1:s1-e1, ctg2:s2-e2, ...`) whenever an MGC spans multiple
  contigs.
- Default `--threads 0` on HPC schedulers (SLURM/PBS cgroups) now resolves the effective
  CPU affinity count via `os.sched_getaffinity(0)` and passes that concrete integer to both
  PyHMMER and Clinker (`-j`), preventing Clinker's `multiprocessing.Pool` from spawning
  workers for all physical host CPUs.

### Changed
- Terminal output at the end of large batch runs (> 10 genomes with hits or > 10 Clinker
  HTML batches) folds per-file path listings into concise summary counts while still
  recording every file path in `magcluster.log`.

## [1.0.1] - 2026-10-21

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
- `magcluster db --help` printed a blank usage line (and `run --help` omitted
  `[OPTIONS]`) with recent typer releases (e.g. 0.27); both usage lines are complete again.
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
- Note on platforms: on Apple silicon (arm64), HMMER's SIMD code can change the last printed
  digit of an E-value compared with Linux x86_64 (e.g. `6.43e-194` vs `6.44e-194`). Hits,
  coordinates, scores and clusters are the same. The tests compare the 1.0.0 reference
  outputs byte for byte on Linux x86_64 and with a float tolerance elsewhere.
- Tests check CLI output independently of terminal colour and width, so they also pass on
  GitHub Actions, where Rich/typer force colour output.
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
