"""Command-line interface for MagCluster v1.0 using Typer and Rich."""

import logging
import os
import platform
import shlex
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import typer
from rich.markup import escape
from rich.table import Table
from . import __version__
from .core.cluster import ClusterEngine
from .core.parser import (
    InputFormatError,
    genome_name_from_path,
    is_genome_file,
    parse_input,
    sanitize_genome_name,
)
from .core.predict import GenePredictor
from .core.scoring import ClusterScorer
from .core.search import DEFAULT_SCORE_CUTOFF, MagnetosomeSearcher
from .data import get_hmm_path, get_profiles_summary_path, load_rules
from .utils.io import export_genome_magnetosome_fasta, write_genome_cluster_tables
from .utils.logger import attach_log_file, console, detach_log_file, detail, log, say
from .viz.clinker_export import (
    ClinkerError,
    export_cluster_to_genbank,
    parse_clinker_opts,
    run_clinker_visualization,
)

HELP_TEXT = f"MagCluster (v{__version__}): Magnetosome Gene Cluster (MGC) identification & visualization."

EPILOG_TEXT = """
Examples:
  # Run MGC analysis on a single GenBank/GBFF genome:
  $ magcluster run genome.gbff -o results/

  # Run MGC analysis on FASTA contigs with 8 CPU threads:
  $ magcluster run contigs.fasta -o results/ -t 8

  # Batch process multiple genomes or a folder of genomes:
  $ magcluster run genomes_dir/ -o results/

  # Inspect built-in Magnetosome HMM database and reference family summary:
  $ magcluster db
"""

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}

app = typer.Typer(
    name="magcluster",
    help=HELP_TEXT,
    epilog=EPILOG_TEXT,
    options_metavar="",
    add_completion=False,
    no_args_is_help=True,
    context_settings=CONTEXT_SETTINGS,
)

def version_callback(value: bool):
    if value:
        console.print(f"[bold cyan]MagCluster[/bold cyan] version [bold green]{__version__}[/bold green]")
        raise typer.Exit()

def show_clinker_help():
    """Display Clinker options guide for use with --clinker-opts."""
    console.print("[bold cyan]=== Clinker Options Guide (for use with --clinker-opts) ===[/bold cyan]\n")

    table = Table(title="Recommended Clinker Options for MagCluster", show_header=True)
    table.add_column("Option", style="bold green")
    table.add_column("Description", style="white")
    table.add_column("Example Usage in MagCluster", style="cyan")

    table.add_row(
        "-na, --no_align",
        "Do not align clusters (only plot gene tracks, fast & clean).",
        '--clinker-opts "-na"'
    )
    table.add_row(
        "-ufo, --use_file_order",
        "Display clusters in the exact order of input files (no tree reordering).",
        '--clinker-opts "-ufo"'
    )
    table.add_row(
        "-mo, --matrix_out <file>",
        "Save cluster similarity distance matrix to a file.",
        '--clinker-opts "-mo matrix.tsv"'
    )
    table.add_row(
        "-o, --output <file>",
        "Save pairwise protein alignment details to a file.",
        '--clinker-opts "-o alignments.tsv"'
    )
    table.add_row(
        "-dc, --decimals <int>",
        "Number of decimal places in alignment scores \\[default: 2].",
        '--clinker-opts "-dc 4"'
    )
    console.print(table)

def clinker_help_callback(value: bool):
    if value:
        show_clinker_help()
        raise typer.Exit()

@app.callback()
def main_callback(
    version: Optional[bool] = typer.Option(
        None,
        "-v", "--version",
        help="Show MagCluster version and exit.",
        callback=version_callback,
        is_eager=True,
    )
):
    """MagCluster command-line entrypoint."""
    pass

# Output files and folders written by `magcluster run` (removed by --force before a re-run).
OWNED_FILE_PATTERNS = [
    "*_mgc_genes.tsv",
    "*_magnetosome_proteins.faa",
    "MGC_align.html",
    "magcluster.log",
    # obsolete names from pre-release builds
    "clinker_colour_map.csv",
    "clinker_gene_functions.csv",
    "mgc_clusters.tsv",
    "presence_absence.tsv",
    "summary.tsv",
    "mgc_genes.tsv",
    "magnetosome_proteins.faa",
]
OWNED_DIRS = {"clusters_gbk": ["*.gbk"], "figures": ["*.svg"]}

LOG_DEPENDENCIES = [
    "biopython", "pyrodigal", "pyhmmer", "pyyaml", "typer", "click", "rich",
    "pygenomeviz", "matplotlib", "clinker", "numpy",
]


def _fail(msg: str, code: int = 1):
    console.print(f"[bold red]Error:[/bold red] {msg}")
    raise typer.Exit(code=code)


def _available_cpus() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def _clean_owned_outputs(outdir: Path) -> List[Path]:
    """Delete outputs of a previous `magcluster run` in *outdir* (used by --force)."""
    removed: List[Path] = []
    for pattern in OWNED_FILE_PATTERNS:
        for p in outdir.glob(pattern):
            if p.is_file() or p.is_symlink():
                p.unlink()
                removed.append(p)
    for sub, patterns in OWNED_DIRS.items():
        d = outdir / sub
        if not d.is_dir():
            continue
        for pattern in patterns:
            for p in d.glob(pattern):
                if p.is_file() or p.is_symlink():
                    p.unlink()
                    removed.append(p)
        try:
            d.rmdir()  # only succeeds if nothing else (user files) is left inside
        except OSError:
            pass
    return removed


def _dependency_versions() -> Dict[str, str]:
    from importlib import metadata
    out = {}
    for name in LOG_DEPENDENCIES:
        try:
            out[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            out[name] = "not installed"
    return out


def _collect_inputs(inputs: List[Path]) -> Tuple[List[Path], List[Tuple[Path, str]], List[Path]]:
    """Expand directories, keep genome files, drop duplicates.

    Returns (files, skipped[(path, reason)], duplicates).
    """
    candidates: List[Path] = []
    skipped: List[Tuple[Path, str]] = []
    for inp in inputs:
        if inp.is_dir():
            for p in sorted(inp.iterdir(), key=lambda x: x.name):
                if p.name.startswith("."):
                    continue
                if not p.is_file():
                    continue
                if is_genome_file(p):
                    candidates.append(p)
                else:
                    skipped.append((p, "not a recognised genome file extension"))
        else:
            candidates.append(inp)

    files: List[Path] = []
    duplicates: List[Path] = []
    seen = set()
    for p in candidates:
        key = p.resolve()
        if key in seen:
            duplicates.append(p)
            continue
        seen.add(key)
        files.append(p)
    return files, skipped, duplicates


def _assign_genome_names(files: List[Path]) -> Tuple[List[str], List[str]]:
    """Unique, file-system safe genome IDs for *files*; returns (names, notices)."""
    names: List[str] = []
    notices: List[str] = []
    used: Dict[str, Path] = {}
    for p in files:
        raw = genome_name_from_path(p)
        safe = sanitize_genome_name(raw)
        if safe != raw:
            notices.append(
                f"Genome name '{raw}' contains spaces or special characters; using '{safe}' "
                f"(input: {p})."
            )
        name = safe
        n = 2
        while name.lower() in used:  # case-insensitive: output files must not collide on macOS/Windows
            name = f"{safe}_{n}"
            n += 1
        if name != safe:
            notices.append(
                f"Genome name '{safe}' is already used by {used[safe.lower()]}; "
                f"{p} will be reported as '{name}'."
            )
        used[name.lower()] = p
        names.append(name)
    return names, notices


@app.command(
    name="run",
    context_settings=CONTEXT_SETTINGS,
    epilog="""
Examples:
  # Analyze a single GenBank/GBFF genome:
  $ magcluster run genome.gbff -o results/

  # Analyze FASTA contigs with 8 CPU threads:
  $ magcluster run contigs.fasta -o results/ -t 8

  # Batch process genomes in a directory:
  $ magcluster run genomes_dir/ -o results/

  # Pass custom options to Clinker (e.g. skip pairwise alignments or keep file order):
  $ magcluster run *.gbff -o results/ --clinker-opts "-na -ufo"

  # Force overwrite when output directory already exists:
  $ magcluster run *.gbff -o results/ --force

  # View available Clinker visualization options:
  $ magcluster run --clinker-help
"""
)
def run_pipeline(
    inputs: List[Path] = typer.Argument(
        ...,
        help="Input genome files (FASTA, GenBank/GBFF, FAA; optionally .gz) or directories.",
        exists=True
    ),
    outdir: Path = typer.Option(
        Path("magcluster_results"),
        "-o", "--outdir",
        help="Output directory."
    ),
    threads: int = typer.Option(
        0,
        "-t", "--threads",
        min=0,
        help="Number of CPU threads (0 = all available cores)."
    ),
    evalue: float = typer.Option(
        1e-5,
        "-e", "--evalue",
        min=0.0,
        help="HMM search E-value cutoff."
    ),
    score: Optional[float] = typer.Option(
        None,
        "-s", "--score",
        show_default=False,
        help="Minimum HMM bit score, applied to all profiles. If omitted, the per-profile "
             "trusted cutoffs from rules.yaml are used (25.0 for all built-in profiles).  "
             "\\[default: per-profile cutoff (25.0)]"
    ),
    max_gap: int = typer.Option(
        5,
        "--max-gap",
        min=0,
        help="Maximum intergenic non-MGC genes allowed between cluster hits."
    ),
    threshold: int = typer.Option(
        3,
        "-th", "--threshold", "--min-families",
        min=1,
        help="Minimum number of distinct magnetosome families (any non-Feo family) required to define an MGC."
    ),
    flank: int = typer.Option(
        3,
        "-un", "--unmag-num", "--flank",
        min=0,
        help="Number of flanking non-MGC genes to reserve on each side of cluster."
    ),
    meta: Optional[bool] = typer.Option(
        None,
        "--meta/--single",
        help="Force metagenome mode vs single-genome mode in Pyrodigal (default: automatic)."
    ),
    identity: float = typer.Option(
        0.3,
        "-i", "--identity",
        min=0.0, max=1.0,
        help="Minimum alignment sequence identity for Clinker (0-1)."
    ),
    clinker_opts: Optional[str] = typer.Option(
        None,
        "--clinker-opts",
        help="Additional options passed directly to clinker (e.g. '--clinker-opts \"-na -ufo\"'). "
             "Relative output paths (e.g. -mo matrix.tsv) are written inside the output directory. "
             "Tip: use '--clinker-help' to see available clinker options."
    ),
    clinker_help: Optional[bool] = typer.Option(
        None,
        "--clinker-help",
        help="Display help and available options for Clinker, then exit.",
        callback=clinker_help_callback,
        is_eager=True,
    ),
    static_plot: bool = typer.Option(
        False,
        "--static-plot",
        help="Also render static vector figures (SVG) with pyGenomeViz into <outdir>/figures/."
    ),
    no_plot: bool = typer.Option(
        False,
        "--no-plot",
        help="Skip Clinker alignment and visualization."
    ),
    force: bool = typer.Option(
        False,
        "-f", "--force",
        help="Overwrite if the output directory already exists and is not empty "
             "(outputs of a previous run there are removed first)."
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        "--log",
        dir_okay=False,
        help="Write the run log to this file instead of <outdir>/magcluster.log."
    ),
):
    """Run full pipeline: Predict -> Search -> Cluster -> Score Mineral -> Visualize with Clinker."""
    start_time = time.time()

    # --- Validate options that Typer cannot check by itself ------------------
    try:
        parse_clinker_opts(clinker_opts)
    except ValueError as e:
        raise typer.BadParameter(str(e), param_hint="'--clinker-opts'")

    console.print(f"[bold cyan]=== MagCluster v{__version__} Execution ===[/bold cyan]")
    console.print(f"[dim]Output directory: {escape(str(outdir.resolve()))}[/dim]\n")

    # --- Collect inputs ------------------------------------------------------
    expanded_files, skipped_files, duplicate_files = _collect_inputs(inputs)

    if not expanded_files:
        for p, reason in skipped_files:
            console.print(f"[yellow]! Skipped {escape(str(p))}: {reason}.[/yellow]")
        console.print("[bold red]No input files found![/bold red]")
        raise typer.Exit(code=1)

    # --- Output directory ----------------------------------------------------
    if outdir.exists() and not outdir.is_dir():
        _fail(f"Output path '[bold]{escape(str(outdir))}[/bold]' exists and is not a directory. "
              f"Please choose a directory with -o/--outdir.")
    removed: List[Path] = []
    if outdir.exists() and any(p for p in outdir.iterdir() if not p.name.startswith(".")):
        if not force:
            console.print(
                f"[bold red]Error:[/bold red] Output directory '[bold]{escape(str(outdir))}[/bold]' already exists and is not empty.\n"
                f"Please specify a different directory or use [bold cyan]--force[/bold cyan] (or [bold cyan]-f[/bold cyan]) to overwrite."
            )
            raise typer.Exit(code=1)
        removed = _clean_owned_outputs(outdir)

    console.print(f"Discovered [bold]{len(expanded_files)}[/bold] input file(s).")
    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        _fail(f"cannot create output directory '{escape(str(outdir))}': {e.strerror or e}")

    log_file = log_path if log_path is not None else outdir / "magcluster.log"
    try:
        log_handler = attach_log_file(log_file)
    except OSError as e:
        _fail(f"cannot write log file '{escape(str(log_file))}': {e.strerror or e}")

    try:
        exit_code = _run_genomes(
            start_time=start_time, inputs=inputs, expanded_files=expanded_files,
            skipped_files=skipped_files, duplicate_files=duplicate_files, removed=removed,
            outdir=outdir, log_file=log_file, threads=threads, evalue=evalue, score=score,
            max_gap=max_gap, threshold=threshold, flank=flank, meta=meta, identity=identity,
            clinker_opts=clinker_opts, static_plot=static_plot, no_plot=no_plot, force=force,
        )
    except BaseException as e:
        if not isinstance(e, (typer.Exit, SystemExit)):
            detail("Run aborted:\n" + traceback.format_exc(), logging.ERROR)
        raise
    finally:
        detach_log_file(log_handler)
    if exit_code:
        raise typer.Exit(code=exit_code)


def _run_genomes(
    *, start_time, inputs, expanded_files, skipped_files, duplicate_files, removed,
    outdir, log_file, threads, evalue, score, max_gap, threshold, flank, meta,
    identity, clinker_opts, static_plot, no_plot, force,
) -> int:
    """Body of `magcluster run` after the output directory and log are ready."""
    # --- Run-log header (file only) ------------------------------------------
    detail(f"MagCluster v{__version__} run log")
    detail(f"Started: {datetime.now().astimezone().isoformat(timespec='seconds')}")
    detail(f"Command: {Path(sys.argv[0]).name} {shlex.join(sys.argv[1:])}")
    detail(f"Working directory: {Path.cwd()}")
    detail(f"Python: {platform.python_version()} ({sys.executable}) on {platform.platform()}")
    detail("Dependencies: " + ", ".join(f"{k} {v}" for k, v in _dependency_versions().items()))

    available = _available_cpus()
    eff_threads = threads
    if threads > available:
        eff_threads = available
    detail("Effective parameters:")
    params = [
        ("inputs", ", ".join(str(p) for p in inputs)),
        ("outdir", str(outdir.resolve())),
        ("log", str(Path(log_file).resolve())),
        ("threads", f"{eff_threads} (0 = all {available} available cores)" if eff_threads == 0 else str(eff_threads)),
        ("evalue", f"{evalue:g}"),
        ("score", f"{score:g} (explicit, applied to all profiles)" if score is not None
                  else f"per-profile trusted cutoff from rules.yaml (default {DEFAULT_SCORE_CUTOFF:g})"),
        ("max_gap", str(max_gap)),
        ("max_distance_bp", "25000 (fixed)"),
        ("threshold (min distinct non-Feo families)", str(threshold)),
        ("flank (unmag-num)", str(flank)),
        ("gene prediction mode", "automatic" if meta is None else ("meta" if meta else "single")),
        ("identity (clinker)", f"{identity:g}"),
        ("clinker_opts", clinker_opts or "-"),
        ("static_plot", str(static_plot)),
        ("no_plot", str(no_plot)),
        ("force", str(force)),
    ]
    for k, v in params:
        detail(f"  {k}: {v}")
    if removed:
        detail(f"--force: removed {len(removed)} output file(s) of a previous run: "
               + ", ".join(str(p.relative_to(outdir)) for p in removed))

    for p, reason in skipped_files:
        say(f"[yellow]! Skipped {escape(str(p))}: {reason}.[/yellow]", logging.WARNING)
    for p in duplicate_files:
        say(f"[yellow]! Skipped duplicate input {escape(str(p))} (already listed).[/yellow]", logging.WARNING)
    if threads > available:
        say(f"[yellow]! --threads {threads} exceeds the {available} available CPU cores; using {available}.[/yellow]",
            logging.WARNING)

    genome_names, name_notices = _assign_genome_names(expanded_files)
    for msg in name_notices:
        say(f"[yellow]! {escape(msg)}[/yellow]", logging.WARNING)
    detail("Inputs (genome ID <- file):")
    for gname, p in zip(genome_names, expanded_files):
        detail(f"  {gname} <- {p.resolve()}")

    # Initialize engines
    predictor = GenePredictor(meta=meta)
    searcher = MagnetosomeSearcher(evalue_cutoff=evalue, score_cutoff=score, cpus=eff_threads)
    clusterer = ClusterEngine(
        max_intergenic_genes=max_gap,
        min_distinct_families=threshold,
        flank_genes=flank
    )

    all_clusters = []
    all_scores = {}
    hits_by_genome: Dict[str, Dict] = {}
    gbk_files: List[Path] = []
    written_tsvs: List[Path] = []
    written_faas: List[Path] = []
    failed: List[Tuple[Path, str]] = []
    summary_rows: List[str] = []
    gbk_dir = outdir / "clusters_gbk"

    n_files = len(expanded_files)
    for idx, (f_path, genome_name) in enumerate(zip(expanded_files, genome_names), 1):
        say(f"[{idx}/{n_files}] Processing [bold green]{escape(f_path.name)}[/bold green]...")
        t0 = time.time()
        try:
            status = _process_one_genome(
                f_path, genome_name, predictor, searcher, clusterer, threshold, meta, outdir,
                gbk_dir, all_clusters, all_scores, hits_by_genome, gbk_files, written_tsvs, written_faas,
            )
        except InputFormatError as e:
            status = f"FAILED: {e}"
            failed.append((f_path, str(e)))
            say(f"  [bold red]Error:[/bold red] cannot read {escape(f_path.name)}: {escape(str(e))}. Skipped.",
                logging.ERROR)
        except Exception as e:  # keep the batch going; details go to the log
            status = f"FAILED: {type(e).__name__}: {e}"
            failed.append((f_path, f"{type(e).__name__}: {e}"))
            say(f"  [bold red]Error:[/bold red] failed to process {escape(f_path.name)} "
                f"({escape(type(e).__name__)}: {escape(str(e))}). Skipped; see the log for details.",
                logging.ERROR)
            detail(traceback.format_exc(), logging.ERROR)
        elapsed_g = time.time() - t0
        detail(f"  {genome_name}: {status} [{elapsed_g:.2f}s]")
        summary_rows.append(f"{genome_name}\t{f_path}\t{elapsed_g:.2f}\t{status}")

    # 5. Export results (per-genome tables and FASTA were written right after each genome)
    say("\n[bold cyan]Writing per-genome reports and magnetosome sequences...[/bold cyan]")

    clinker_ok = False
    clinker_failed = False
    plot_failed = False
    figure_files: List[Path] = []
    # 6. Clinker visualization
    if all_clusters:
        say("[dim]-> Preparing sliced GenBank files for Clinker visualization...[/dim]")

        if not no_plot:
            clinker_html = outdir / "MGC_align.html"
            say(f"[dim]-> Running Clinker to generate interactive comparison: {clinker_html.name}...[/dim]")
            hits_by_cluster = {c.cluster_id: hits_by_genome.get(c.genome_name, {}) for c in all_clusters}
            try:
                clinker_ok = run_clinker_visualization(
                    gbk_files=gbk_files,
                    plot_html=clinker_html,
                    hits={},
                    jobs=eff_threads,
                    identity=identity,
                    clinker_opts=clinker_opts,
                    hits_by_cluster=hits_by_cluster,
                    cwd=outdir,
                    raise_on_error=True,
                )
            except ClinkerError as e:
                clinker_failed = True
                say(f"  [bold red]Error:[/bold red] {escape(str(e))}", logging.ERROR)
            if clinker_ok:
                say(f"  [bold green]✓ Clinker interactive alignment created: {escape(str(clinker_html))}[/bold green]")

        if static_plot:
            say(f"[dim]-> Rendering static cluster maps with pyGenomeViz...[/dim]")
            try:
                from .viz.static_plot import plot_all_clusters
                figure_files = plot_all_clusters(all_clusters, {}, outdir / "figures", fmt="svg",
                                                 hits_by_genome=hits_by_genome)
            except Exception as e:
                plot_failed = True
                say(f"  [bold red]Error:[/bold red] static figure rendering failed "
                    f"({escape(type(e).__name__)}: {escape(str(e))}).", logging.ERROR)
                detail(traceback.format_exc(), logging.ERROR)

    elapsed = time.time() - start_time
    n_ok = n_files - len(failed)
    if failed and n_ok == 0:
        say(f"\n[bold red]✗ No input could be processed ({len(failed)} failed).[/bold red] See the log: {escape(str(log_file))}",
            logging.ERROR)
    elif failed or clinker_failed or plot_failed:
        problems = []
        if failed:
            problems.append(f"{len(failed)} of {n_files} input(s) failed")
        if clinker_failed:
            problems.append("Clinker failed")
        if plot_failed:
            problems.append("static plotting failed")
        say(f"\n[bold yellow]! Finished with errors in {elapsed:.1f}s ({'; '.join(problems)}).[/bold yellow] "
            f"Results saved to: [bold]{escape(str(outdir.resolve()))}[/bold]", logging.WARNING)
    else:
        say(f"\n[bold green]✓ Done in {elapsed:.1f}s![/bold green] Results saved to: [bold]{escape(str(outdir.resolve()))}[/bold]")
    for tsv_p in written_tsvs:
        say(f"  • Magnetosome Genes: [underline]{escape(str(tsv_p))}[/underline]")
    for faa_p in written_faas:
        say(f"  • Magnetosome Proteins: [underline]{escape(str(faa_p))}[/underline]")
    if clinker_ok:
        say(f"  • Clinker Visualization: [underline]{escape(str(outdir))}/MGC_align.html[/underline]")
    for fig in figure_files:
        detail(f"  • Static figure: {fig}")
    for p in gbk_files:
        detail(f"  • Cluster GenBank: {p}")
    console.print(f"  • Run Log: [underline]{escape(str(log_file))}[/underline]")

    detail("Per-genome summary (genome_id, input, seconds, result):")
    for row in summary_rows:
        detail("  " + row)
    for p, reason in failed:
        detail(f"Failed input: {p}: {reason}", logging.ERROR)
    code = 1 if (failed or clinker_failed or plot_failed) else 0
    detail(f"Finished: {datetime.now().astimezone().isoformat(timespec='seconds')} "
           f"(elapsed {elapsed:.1f}s, exit code {code})")
    return code


def _process_one_genome(
    f_path, genome_name, predictor, searcher, clusterer, threshold, meta, outdir, gbk_dir,
    all_clusters, all_scores, hits_by_genome, gbk_files, written_tsvs, written_faas,
) -> str:
    """Parse, predict, search and cluster one input; write its outputs. Returns a status line."""
    # 1. Parse input
    parsed = parse_input(f_path, genome_name=genome_name)
    n_bp = sum(len(r.seq) for r in parsed.contigs.values())
    detail(f"  {genome_name}: format={parsed.file_type}, records={len(parsed.contigs)}, bp={n_bp}")

    # 2. Predict genes if needed
    if parsed.needs_prediction:
        say(f"  [dim]-> Predicting genes with Pyrodigal...[/dim]")
        genes = predictor.predict_from_records(list(parsed.contigs.values()), genome_name=genome_name)
        detail(f"  {genome_name}: Pyrodigal mode={predictor.last_mode}, genes={len(genes)}")
        if predictor.last_fallback_reason and meta is False:
            say(f"  [yellow]! --single: {escape(predictor.last_fallback_reason)}.[/yellow]", logging.WARNING)
        elif predictor.last_fallback_reason:
            detail(f"  {genome_name}: {predictor.last_fallback_reason}")
    else:
        genes = parsed.genes
        say(f"  [dim]-> Found {len(genes)} pre-annotated CDS.[/dim]")
        if meta is not None:
            flag = "--meta" if meta else "--single"
            say(f"  [yellow]! {flag} ignored for {escape(f_path.name)}: it already contains CDS/protein "
                f"annotations, so no gene prediction is run.[/yellow]", logging.WARNING)

    genes_lookup = {}
    for g in genes:
        g.genome_name = genome_name
        genes_lookup[g.gene_id] = g

    if not genes:
        say(f"  [yellow]No genes found in {escape(f_path.name)}, skipping.[/yellow]", logging.WARNING)
        return "no genes found"

    # 3. Search HMMs
    say(f"  [dim]-> Searching {len(genes)} proteins against magnetosome HMM library...[/dim]")
    hits = searcher.search_genes(genes)
    hits_by_genome[genome_name] = hits
    say(f"  [dim]-> Identified {len(hits)} significant magnetosome gene hits.[/dim]")

    if not hits:
        return f"{len(genes)} genes, 0 magnetosome hits"

    # 4. Cluster
    clusters = []
    if not parsed.has_coordinates:
        say(f"  [yellow]! Input protein FASTA lacks genomic coordinates (contig/start/end).[/yellow]", logging.WARNING)
        say(f"    [dim]Identified {len(hits)} magnetosome protein(s); syntenic cluster detection and Clinker plotting skipped.[/dim]")
    else:
        contig_lengths = {cid: len(rec.seq) for cid, rec in parsed.contigs.items()}
        clusters = clusterer.cluster_genome(
            genome_name, genes, hits,
            contig_lengths=contig_lengths, circular_contigs=parsed.circular_contigs,
        )
        if clusters:
            say(f"  [bold cyan]-> Detected {len(clusters)} MGC cluster(s)![/bold cyan]")
            for c in clusters:
                score_obj = ClusterScorer.score_cluster(c)
                all_clusters.append(c)
                all_scores[c.cluster_id] = score_obj
                sub_count = len(c.subclusters) if getattr(c, "subclusters", None) else 1
                sub_info = f" ({sub_count} subclusters, {len(c.all_genes)} genes)" if sub_count > 1 else f" ({len(c.all_genes)} genes)"
                origin = " (spans origin)" if c.wraps_origin else ""
                say(
                    f"     • [bold]{escape(c.cluster_id)}[/bold]{sub_info}: {escape(c.contig_id)}:{c.start_bp}-{c.end_bp}{origin} "
                    f"| Families: {score_obj.n_mag_families} "
                    f"| Mineral: [bold]{score_obj.mineral_type}[/bold]"
                )
                for sub in c.subclusters:
                    fams = sorted(sub.distinct_magnetosome_families)
                    detail(f"       {sub.subcluster_id}: {sub.contig_id}:{sub.start_bp}-{sub.end_bp} "
                           f"({len(sub.all_genes)} genes; {len(fams)} families: {', '.join(fams)})"
                           f"{' [contig edge]' if sub.on_edge else ''}")
                detail(f"       families: {', '.join(score_obj.magnetosome_families_found)}; "
                       f"flags: {', '.join(score_obj.flags) or '-'}")
        else:
            say(f"  [dim]-> No syntenic MGC clusters met the threshold criteria (>= {threshold} distinct magnetosome families).[/dim]")

    # Write this genome's outputs immediately, so later failures cannot lose them
    written_tsvs.append(write_genome_cluster_tables(genome_name, clusters, hits, outdir, genes_lookup=genes_lookup))
    written_faas.append(export_genome_magnetosome_fasta(genome_name, clusters, hits, outdir, genes_lookup=genes_lookup))
    for c in clusters:
        gbk_file = gbk_dir / f"{c.cluster_id}.gbk"
        export_cluster_to_genbank(c, parsed.contigs, hits, gbk_file)
        gbk_files.append(gbk_file)

    n_clustered = len({g.gene_id for c in clusters for g in c.all_genes if g.gene_id in hits})
    if clusters:
        c = clusters[0]
        sc = all_scores[c.cluster_id]
        return (f"{len(genes)} genes, {len(hits)} hits ({n_clustered} in MGC); MGC {c.cluster_id} "
                f"{c.contig_id}:{c.start_bp}-{c.end_bp}{' (spans origin)' if c.wraps_origin else ''}, "
                f"{len(c.subclusters)} subcluster(s), {sc.n_mag_families} families, {sc.mineral_type}")
    if not parsed.has_coordinates:
        return f"{len(genes)} proteins, {len(hits)} hits (no coordinates, clustering skipped)"
    return f"{len(genes)} genes, {len(hits)} hits, no MGC"


@app.command(name="db", context_settings=CONTEXT_SETTINGS)
def db_info():
    """Display information about the built-in Magnetosome HMM database and rules."""
    import csv

    summary_p = get_profiles_summary_path()
    rules = load_rules()
    user_db_path = get_hmm_path().resolve()

    console.print("[bold cyan]=== MagCluster Reference Database Info ===[/bold cyan]")
    console.print(f"Rules version: [bold green]{rules.get('version', 'unknown')}[/bold green]")
    console.print(f"HMM database path: [dim]{user_db_path}[/dim]")
    console.print(f"Total HMM profiles: [bold]{len(rules.get('profiles', {}))}[/bold]")
    console.print(f"Core defining families ({len(rules.get('core_families', []))}): [bold magenta]{', '.join(rules.get('core_families', []))}[/bold magenta]\n")

    if summary_p.exists():
        table = Table(title="Reference Magnetosome HMM Profiles (92)", show_header=True)
        table.add_column("Profile", style="bold cyan")
        table.add_column("Protein Family", style="bold green")
        table.add_column("Ref Seqs", justify="right", style="cyan")
        table.add_column("Median length (aa)", justify="right", style="magenta")

        with open(summary_p, "r", encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))

        # Sort descending by Ref Seqs count, then alphabetically by Profile name
        sorted_rows = sorted(rows, key=lambda r: (-int(r["n_seqs"]), r["profile"]))

        for r in sorted_rows:
            table.add_row(
                r["profile"],
                r["parent_family"],
                r["n_seqs"],
                r["median_length"],
            )

        console.print(table)

if __name__ == "__main__":
    app()
