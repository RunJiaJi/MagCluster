"""Command-line interface for MagCluster v1.0 using Typer and Rich."""

import sys
import time
from pathlib import Path
from typing import List, Optional
import typer
from rich.table import Table
from . import __version__
from .core.cluster import ClusterEngine
from .core.parser import parse_input
from .core.predict import GenePredictor
from .core.scoring import ClusterScorer
from .core.search import MagnetosomeSearcher
from .data import get_hmm_path, get_profiles_summary_path, load_rules
from .utils.io import export_all_genome_fastas, write_all_genome_tables
from .utils.logger import console, log
from .viz.clinker_export import export_all_clusters_for_clinker, run_clinker_visualization
from .viz.static_plot import plot_all_clusters

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

app = typer.Typer(
    name="magcluster",
    help=HELP_TEXT,
    epilog=EPILOG_TEXT,
    options_metavar="",
    add_completion=False,
    no_args_is_help=True
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
        "Number of decimal places in alignment scores [default: 2].",
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

@app.command(
    name="run",
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
        help="Input genome files (FASTA, GenBank/GBFF, FAA) or directories.",
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
        help="Number of CPU threads (0 = all available cores)."
    ),
    evalue: float = typer.Option(
        1e-5,
        "-e", "--evalue",
        help="HMM search E-value cutoff."
    ),
    score: float = typer.Option(
        25.0,
        "-s", "--score",
        help="Minimum bit score threshold for HMM hits."
    ),
    max_gap: int = typer.Option(
        5,
        "--max-gap",
        help="Maximum intergenic non-MGC genes allowed between cluster hits."
    ),
    threshold: int = typer.Option(
        3,
        "-th", "--threshold", "--min-families",
        help="Minimum number of distinct magnetosome families required to define an MGC (default: 3)."
    ),
    flank: int = typer.Option(
        3,
        "-un", "--unmag-num", "--flank",
        help="Number of flanking non-MGC genes to reserve on each side of cluster (default: 3)."
    ),
    meta: Optional[bool] = typer.Option(
        None,
        "--meta/--single",
        help="Force metagenome mode vs single-genome mode in Pyrodigal."
    ),
    identity: float = typer.Option(
        0.3,
        "-i", "--identity",
        help="Minimum alignment sequence identity for Clinker (default: 0.3)."
    ),
    clinker_opts: Optional[str] = typer.Option(
        None,
        "--clinker-opts",
        help="Additional options passed directly to clinker (e.g. '--clinker-opts \"-na -ufo\"'). Tip: use '--clinker-help' to see available clinker options."
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
        help="Also render static vector figures (SVG) with pyGenomeViz."
    ),
    no_plot: bool = typer.Option(
        False,
        "--no-plot",
        help="Skip Clinker alignment and visualization."
    ),
    force: bool = typer.Option(
        False,
        "-f", "--force",
        help="Overwrite existing files if the output directory already exists and is not empty."
    ),
):
    """Run full pipeline: Predict -> Search -> Cluster -> Score Mineral -> Visualize with Clinker."""
    start_time = time.time()
    console.print(f"[bold cyan]=== MagCluster v{__version__} Execution ===[/bold cyan]")
    console.print(f"[dim]Output directory: {outdir.resolve()}[/dim]\n")

    # Expand inputs if directory
    expanded_files = []
    for inp in inputs:
        if inp.is_dir():
            for p in inp.iterdir():
                if p.is_file() and not p.name.startswith("."):
                    expanded_files.append(p)
        else:
            expanded_files.append(inp)

    if not expanded_files:
        console.print("[bold red]No input files found![/bold red]")
        raise typer.Exit(code=1)

    # Check if outdir exists and is non-empty
    if outdir.exists() and any(p for p in outdir.iterdir() if not p.name.startswith(".")):
        if not force:
            console.print(
                f"[bold red]Error:[/bold red] Output directory '[bold]{outdir}[/bold]' already exists and is not empty.\n"
                f"Please specify a different directory or use [bold cyan]--force[/bold cyan] (or [bold cyan]-f[/bold cyan]) to overwrite."
            )
            raise typer.Exit(code=1)

    console.print(f"Discovered [bold]{len(expanded_files)}[/bold] input file(s).")
    outdir.mkdir(parents=True, exist_ok=True)

    # Clean up obsolete table files if re-running in existing outdir
    for old_file in [
        "clinker_colour_map.csv",
        "clinker_gene_functions.csv",
        "mgc_clusters.tsv",
        "presence_absence.tsv",
        "summary.tsv",
        "mgc_genes.tsv",
        "magnetosome_proteins.faa",
    ]:
        p = outdir / old_file
        if p.exists():
            p.unlink()

    # Initialize engines
    predictor = GenePredictor(meta=meta)
    searcher = MagnetosomeSearcher(evalue_cutoff=evalue, score_cutoff=score, cpus=threads)
    clusterer = ClusterEngine(
        max_intergenic_genes=max_gap,
        min_distinct_families=threshold,
        flank_genes=flank
    )

    all_clusters = []
    all_scores = {}
    all_hits = {}
    all_genes_lookup = {}
    contigs_lookup = {}
    processed_genomes: List[str] = []

    for idx, f_path in enumerate(expanded_files, 1):
        console.print(f"[{idx}/{len(expanded_files)}] Processing [bold green]{f_path.name}[/bold green]...")
        
        # 1. Parse input
        parsed = parse_input(f_path)
        genome_name = parsed.name
        contigs_lookup[genome_name] = parsed.contigs

        # 2. Predict genes if needed
        if parsed.needs_prediction:
            console.print(f"  [dim]-> Predicting genes with Pyrodigal...[/dim]")
            genes = predictor.predict_from_records(list(parsed.contigs.values()), genome_name=genome_name)
        else:
            genes = parsed.genes
            console.print(f"  [dim]-> Found {len(genes)} pre-annotated CDS.[/dim]")

        for g in genes:
            g.genome_name = genome_name
            all_genes_lookup[g.gene_id] = g

        if not genes:
            console.print(f"  [yellow]No genes found in {f_path.name}, skipping.[/yellow]")
            continue

        # 3. Search HMMs
        console.print(f"  [dim]-> Searching {len(genes)} proteins against magnetosome HMM library...[/dim]")
        hits = searcher.search_genes(genes)
        all_hits.update(hits)
        console.print(f"  [dim]-> Identified {len(hits)} significant magnetosome gene hits.[/dim]")

        if not hits:
            continue

        if genome_name not in processed_genomes:
            processed_genomes.append(genome_name)

        # 4. Cluster
        if not parsed.has_coordinates:
            console.print(f"  [yellow]! Input protein FASTA lacks genomic coordinates (contig/start/end).[/yellow]")
            console.print(f"    [dim]Identified {len(hits)} magnetosome protein(s); syntenic cluster detection and Clinker plotting skipped.[/dim]")
            continue

        clusters = clusterer.cluster_genome(genome_name, genes, hits)
        if clusters:
            console.print(f"  [bold cyan]-> Detected {len(clusters)} MGC cluster(s)![/bold cyan]")
            for c in clusters:
                score_obj = ClusterScorer.score_cluster(c)
                all_clusters.append(c)
                all_scores[c.cluster_id] = score_obj
                sub_count = len(c.subclusters) if getattr(c, "subclusters", None) else 1
                sub_info = f" ({sub_count} subclusters, {len(c.all_genes)} genes)" if sub_count > 1 else f" ({len(c.all_genes)} genes)"
                console.print(
                    f"     • [bold]{c.cluster_id}[/bold]{sub_info}: {c.contig_id}:{c.start_bp}-{c.end_bp} "
                    f"| Families: {score_obj.n_mag_families} "
                    f"| Mineral: [bold]{score_obj.mineral_type}[/bold]"
                )
        else:
            console.print(f"  [dim]-> No syntenic MGC clusters met the threshold criteria (>= {threshold} distinct magnetosome families).[/dim]")

    # 5. Export results
    console.print("\n[bold cyan]Writing per-genome reports and magnetosome sequences...[/bold cyan]")
    written_tsvs = write_all_genome_tables(processed_genomes, all_clusters, all_scores, all_hits, outdir, genes_lookup=all_genes_lookup)
    written_faas = export_all_genome_fastas(processed_genomes, all_clusters, all_hits, outdir, genes_lookup=all_genes_lookup)

    # 6. Clinker visualization
    if all_clusters:
        console.print("[dim]-> Preparing sliced GenBank files for Clinker visualization...[/dim]")
        gbk_files = export_all_clusters_for_clinker(all_clusters, contigs_lookup, all_hits, outdir)

        if not no_plot:
            clinker_html = outdir / "MGC_align.html"
            console.print(f"[dim]-> Running Clinker to generate interactive comparison: {clinker_html.name}...[/dim]")
            success = run_clinker_visualization(
                gbk_files=gbk_files,
                plot_html=clinker_html,
                hits=all_hits,
                jobs=threads,
                identity=identity,
                clinker_opts=clinker_opts
            )
            if success:
                console.print(f"  [bold green]✓ Clinker interactive alignment created: {clinker_html}[/bold green]")

        if static_plot:
            console.print(f"[dim]-> Rendering static cluster maps with pyGenomeViz...[/dim]")
            plot_all_clusters(all_clusters, all_hits, outdir / "figures", fmt="svg")

    elapsed = time.time() - start_time
    console.print(f"\n[bold green]✓ Done in {elapsed:.1f}s![/bold green] Results saved to: [bold]{outdir.resolve()}[/bold]")
    if written_tsvs:
        for tsv_p in written_tsvs:
            console.print(f"  • Magnetosome Genes: [underline]{tsv_p}[/underline]")
    if written_faas:
        for faa_p in written_faas:
            console.print(f"  • Magnetosome Proteins: [underline]{faa_p}[/underline]")
    if all_clusters and not no_plot:
        console.print(f"  • Clinker Visualization: [underline]{outdir}/MGC_align.html[/underline]")

@app.command(name="db")
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
