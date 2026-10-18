"""Publication-quality cluster visualization using pyGenomeViz (matplotlib-based)."""

from pathlib import Path
from typing import Dict, List, Optional
from pygenomeviz import GenomeViz
from ..core.cluster import MGCCluster
from ..core.search import MagnetosomeHit

def plot_mgc_cluster(
    cluster: MGCCluster,
    hits: Dict[str, MagnetosomeHit],
    out_path: Path,
    fig_width: int = 14,
    feature_style: str = "arrow"
):
    """Render a single MGC cluster to SVG, PNG, or PDF."""
    gv = GenomeViz(
        fig_track_height=1.0,
        feature_track_ratio=0.8,
        link_track_ratio=0.8,
    )

    track_name = f"{cluster.cluster_id} ({cluster.contig_id})"
    track_size = cluster.span_bp
    track = gv.add_feature_track(track_name, track_size)

    # Offset to make plot relative to cluster start_bp
    offset = cluster.start_bp

    for g in cluster.all_genes:
        g_start = max(0, g.start - offset)
        g_end = min(track_size, g.end - offset)
        strand = g.strand

        # Determine color and label
        if g.gene_id in hits:
            hit = hits[g.gene_id]
            color = hit.color
            label = hit.profile[:-5] if hit.profile.endswith("_full") else hit.profile
        else:
            color = "#E0E0E0"  # gray for intervening non-MGC genes
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
    fmt: str = "svg"
):
    """Plot all clusters in the list."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for c in clusters:
        f_name = out_dir / f"{c.cluster_id}.{fmt}"
        plot_mgc_cluster(c, hits, f_name)
