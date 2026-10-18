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

app = typer.Typer(
    name="magcluster",
    help="MagCluster: Fast, Prokka-free Magnetosome Gene Cluster (MGC) identification & visualization.",
    add_completion=False,
    no_args_is_help=True
)

@app.command(name="version")
def version_cmd():
    """Display MagCluster version."""
    console.print(f"[bold cyan]MagCluster[/bold cyan] version [bold green]{__version__}[/bold green]")

@app.command(name="run")
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
        console.print("[dim]-> Preparing sliced GenBank files and Clinker gene function maps...[/dim]")
        gbk_files, gf_file, cm_file = export_all_clusters_for_clinker(all_clusters, contigs_lookup, all_hits, outdir)

        if not no_plot:
            clinker_html = outdir / "MGC_align.html"
            console.print(f"[dim]-> Running Clinker to generate interactive comparison: {clinker_html.name}...[/dim]")
            success = run_clinker_visualization(
                gbk_files=gbk_files,
                gene_functions_csv=gf_file,
                colour_map_csv=cm_file,
                plot_html=clinker_html,
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
        console.print(f"  • Summary: [underline]{outdir}/summary.tsv[/underline]")
        console.print(f"  • Cluster Table: [underline]{outdir}/mgc_clusters.tsv[/underline]")
        console.print(f"  • Combined Proteins: [underline]{outdir}/magnetosome_proteins.faa[/underline]")
        if not no_plot:
            console.print(f"  • Clinker Figure: [underline]{outdir}/MGC_align.html[/underline]")

@app.command(name="db")
def db_info():
    """Display information about the built-in Magnetosome HMM database and rules."""
    summary_p = get_profiles_summary_path()
    rules = load_rules()

    console.print("[bold cyan]=== MagCluster Reference Database Info ===[/bold cyan]")
    console.print(f"Rules version: [bold green]{rules.get('version', 'unknown')}[/bold green]")
    console.print(f"HMM database path: [dim]{get_hmm_path()}[/dim]")
    console.print(f"Total HMM profiles: [bold]{len(rules.get('profiles', {}))}[/bold]")
    console.print(f"Core defining families ({len(rules.get('core_families', []))}): [bold magenta]{', '.join(rules.get('core_families', []))}[/bold magenta]\n")

    if summary_p.exists():
        table = Table(title="Top 20 Reference Profiles Summary", show_header=True)
        table.add_column("Profile", style="cyan")
        table.add_column("Parent Family", style="green")
        table.add_column("Role", style="magenta")
        table.add_column("Lineage", style="blue")
        table.add_column("Ref Seqs", justify="right")
        table.add_column("Median AA", justify="right")

        with open(summary_p, "r", encoding="utf-8") as f:
            lines = f.readlines()
            for line in lines[1:21]:
                parts = line.strip().split("\t")
                table.add_row(parts[0], parts[1], parts[2], parts[3], parts[4], parts[5])
        console.print(table)

if __name__ == "__main__":
    app()
