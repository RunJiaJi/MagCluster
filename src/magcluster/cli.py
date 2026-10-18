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
from .utils.io import export_combined_magnetosome_fasta, write_cluster_tables
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

    console.print(f"Discovered [bold]{len(expanded_files)}[/bold] input file(s).")
    outdir.mkdir(parents=True, exist_ok=True)

    # Clean up obsolete table files if re-running in existing outdir
    for old_file in [
        "clinker_colour_map.csv",
        "clinker_gene_functions.csv",
        "mgc_clusters.tsv",
        "presence_absence.tsv",
        "summary.tsv",
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
    contigs_lookup = {}

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

        # 4. Cluster
        clusters = clusterer.cluster_genome(genome_name, genes, hits)
        if clusters:
            console.print(f"  [bold cyan]-> Detected {len(clusters)} MGC cluster(s)![/bold cyan]")
            for c in clusters:
                score_obj = ClusterScorer.score_cluster(c)
                all_clusters.append(c)
                all_scores[c.cluster_id] = score_obj
                console.print(
                    f"     • [bold]{c.cluster_id}[/bold]: {c.contig_id}:{c.start_bp}-{c.end_bp} "
                    f"({c.span_bp} bp) | Families: {score_obj.n_mag_families} "
                    f"| Mineral: [bold]{score_obj.mineral_type}[/bold]"
                )
        else:
            console.print(f"  [dim]-> No syntenic MGC clusters met the threshold criteria (>= {threshold} distinct magnetosome families).[/dim]")

    # 5. Export results
    console.print("\n[bold cyan]Writing reports and magnetosome sequences...[/bold cyan]")
    write_cluster_tables(all_clusters, all_scores, all_hits, outdir)
    export_combined_magnetosome_fasta(all_clusters, all_hits, outdir, filename="magnetosome_proteins.faa")

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
                identity=identity
            )
            if success:
                console.print(f"  [bold green]✓ Clinker interactive alignment created: {clinker_html}[/bold green]")

        if static_plot:
            console.print(f"[dim]-> Rendering static cluster maps with pyGenomeViz...[/dim]")
            plot_all_clusters(all_clusters, all_hits, outdir / "figures", fmt="svg")

    elapsed = time.time() - start_time
    console.print(f"\n[bold green]✓ Done in {elapsed:.1f}s![/bold green] Results saved to: [bold]{outdir.resolve()}[/bold]")
    if all_clusters:
        console.print(f"  • Magnetosome Genes: [underline]{outdir}/mgc_genes.tsv[/underline]")
        console.print(f"  • Combined Proteins: [underline]{outdir}/magnetosome_proteins.faa[/underline]")
        if not no_plot:
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
