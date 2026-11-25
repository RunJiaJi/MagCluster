"""Publication-quality cluster visualization using pyGenomeViz (matplotlib-based)."""

from pathlib import Path
from typing import Dict, List, Optional
from ..core.cluster import MGCCluster
from ..core.search import MagnetosomeHit
from ..utils.logger import quiet_matplotlib
from .clinker_export import FAMILY_COLOR_PALETTE, get_display_label, get_family_class


def _genome_viz_class():
    """Import pyGenomeViz (and thus matplotlib) lazily, only when a figure is drawn."""
    quiet_matplotlib()
    from pygenomeviz import GenomeViz
    return GenomeViz

def plot_mgc_cluster(
    cluster: MGCCluster,
    hits: Dict[str, MagnetosomeHit],
    out_path: Path,
    fig_width: int = 14,
    feature_style: str = "arrow"
):
    """Render a single MGC cluster to SVG, PNG, or PDF."""
    GenomeViz = _genome_viz_class()
    gv = GenomeViz(
        fig_width=fig_width,
        fig_track_height=1.0,
        feature_track_ratio=0.8,
        link_track_ratio=0.8,
    )

    subclusters = getattr(cluster, "subclusters", None) or []
    slices = subclusters if subclusters else [cluster]

    for i, s in enumerate(slices, 1):
        ctg_label = s.contig_id
        track_name = f"{cluster.cluster_id} loc{i} ({ctg_label})" if len(slices) > 1 else f"{cluster.cluster_id} ({ctg_label})"
        track_size = s.span_bp
        track = gv.add_feature_track(track_name, track_size)

        offset = s.start_bp

        for g in s.all_genes:
            g_start = max(0, g.start - offset)
            g_end = min(track_size, g.end - offset + 1)
            if g_end <= g_start:
                continue
            strand = g.strand

            # Determine unified 8-family color and label matching Clinker HTML
            if g.gene_id in hits:
                hit = hits[g.gene_id]
                f_class = get_family_class(hit.parent_family, hit)
                color = FAMILY_COLOR_PALETTE.get(f_class, FAMILY_COLOR_PALETTE["Other"])
                label = get_display_label(hit)
            else:
                color = FAMILY_COLOR_PALETTE["Other"]
                label = ""

            text_kws = {
                "size": 8,
                "rotation": 45 if len(label) > 4 else 0,
                "vpos": "top"
            }

            track.add_feature(
                start=g_start,
                end=g_end,
                strand=strand,
                plotstyle=feature_style,
                label=label,
                text_kws=text_kws,
                facecolor=color,
                edgecolor="#333333",
                linewidth=0.5
            )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    gv.savefig(str(out_path))

def plot_all_clusters(
    clusters: List[MGCCluster],
    hits: Dict[str, MagnetosomeHit],
    out_dir: Path,
    fmt: str = "svg",
    hits_by_genome: Optional[Dict[str, Dict[str, MagnetosomeHit]]] = None,
) -> List[Path]:
    """Plot all clusters in the list; returns the written figure paths.

    ``hits_by_genome`` (genome name -> hits of that genome) takes precedence over
    the flat ``hits`` dict so that identical gene IDs in different genomes never mix.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for c in clusters:
        f_name = out_dir / f"{c.cluster_id}.{fmt}"
        g_hits = hits_by_genome.get(c.genome_name, {}) if hits_by_genome is not None else hits
        plot_mgc_cluster(c, g_hits, f_name)
        written.append(f_name)
    return written
