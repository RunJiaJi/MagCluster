"""Clinker-compatible GenBank slice exporter and automated visualization runner."""

import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from ..core.cluster import MGCCluster
from ..core.search import MagnetosomeHit
from ..utils.logger import log

def export_cluster_to_genbank(
    cluster: MGCCluster,
    contig_record: Optional[SeqRecord],
    hits: Dict[str, MagnetosomeHit],
    out_path: Path
):
    """Export a sliced GenBank file for a cluster."""
    offset = cluster.start_bp - 1  # 0-based slice start
    end_idx = cluster.end_bp       # 0-based slice end (exclusive)

    if contig_record and len(contig_record.seq) >= end_idx:
        sub_seq = contig_record.seq[offset:end_idx]
    else:
        # Fallback if no nucleotide contig is available: approximate sequence with Ns
        span_len = cluster.span_bp
        sub_seq = Seq("N" * span_len)

    new_rec = SeqRecord(
        sub_seq,
        id=cluster.cluster_id,
        name=cluster.cluster_id[:16],
        description=f"Magnetosome Gene Cluster from {cluster.genome_name} ({cluster.contig_id}:{cluster.start_bp}-{cluster.end_bp})"
    )
    new_rec.annotations["molecule_type"] = "DNA"
    new_rec.annotations["topology"] = "linear"

    features = []
    # Add source feature
    source_feat = SeqFeature(
        FeatureLocation(0, len(sub_seq)),
        type="source",
        qualifiers={"organism": [cluster.genome_name]}
    )
    features.append(source_feat)

    # Add CDS features
    for g in cluster.all_genes:
        feat_start = max(0, g.start - 1 - offset)
        feat_end = min(len(sub_seq), g.end - offset)
        feat_strand = g.strand

        if feat_end <= feat_start:
            continue

        loc = FeatureLocation(feat_start, feat_end, strand=feat_strand)

        qualifiers = {
            "locus_tag": [g.gene_id],
            "translation": [g.translation]
        }

        if g.gene_id in hits:
            h = hits[g.gene_id]
            qualifiers["product"] = [f"Magnetosome protein {h.parent_family}"]
            qualifiers["gene"] = [h.parent_family]
            qualifiers["note"] = [f"Profile={h.profile}; Score={h.score:.1f}; Evalue={h.evalue:.2e}"]
        else:
            qualifiers["product"] = [g.product or "hypothetical protein"]

        feat = SeqFeature(loc, type="CDS", qualifiers=qualifiers)
        features.append(feat)

    new_rec.features = features
    out_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(new_rec, str(out_path), "genbank")

def export_all_clusters_for_clinker(
    clusters: List[MGCCluster],
    contigs_lookup: Dict[str, Dict[str, SeqRecord]],
    hits: Dict[str, MagnetosomeHit],
    out_dir: Path
) -> Tuple[List[Path], Path, Path]:
    """Export sliced GenBank files, clinker -gf gene_functions.csv, and -cm colour_map.csv."""
    gbk_dir = out_dir / "clusters_gbk"
    gbk_dir.mkdir(parents=True, exist_ok=True)

    gf_lines = []
    colour_map = {}
    gbk_files = []

    for c in clusters:
        gbk_file = gbk_dir / f"{c.cluster_id}.gbk"
        ctg_rec = contigs_lookup.get(c.genome_name, {}).get(c.contig_id)
        export_cluster_to_genbank(c, ctg_rec, hits, gbk_file)
        gbk_files.append(gbk_file)

        for g in c.all_genes:
            if g.gene_id in hits:
                h = hits[g.gene_id]
                func = h.parent_family
                color = h.color
            else:
                func = "Other"
                color = "#BDBDBD"

            gf_lines.append(f"{g.gene_id},{func}\n")
            colour_map[func] = color

    # 1. gene_functions.csv (gene,function)
    gf_file = out_dir / "clinker_gene_functions.csv"
    with open(gf_file, "w", encoding="utf-8") as f:
        f.writelines(gf_lines)

    # 2. colour_map.csv (function,colour)
    cm_file = out_dir / "clinker_colour_map.csv"
    with open(cm_file, "w", encoding="utf-8") as f:
        for func, color in colour_map.items():
            f.write(f"{func},{color}\n")

    return gbk_files, gf_file, cm_file

def run_clinker_visualization(
    gbk_files: List[Path],
    gene_functions_csv: Path,
    colour_map_csv: Path,
    plot_html: Path,
    jobs: int = 0,
    identity: float = 0.3
) -> bool:
    """Run Clinker to generate publication-quality interactive HTML comparison figure."""
    if not gbk_files:
        return False

    cmd = [
        "clinker",
        *[str(p) for p in gbk_files],
        "-gf", str(gene_functions_csv),
        "-cm", str(colour_map_csv),
        "-p", str(plot_html),
        "-j", str(jobs),
        "-i", str(identity)
    ]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return True
    except FileNotFoundError:
        log.warning("clinker is not installed on PATH; skipped interactive HTML figure generation.")
        return False
    except subprocess.CalledProcessError as e:
        log.warning(f"clinker execution warning: {e.stderr}")
        return False
