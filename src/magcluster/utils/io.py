"""Tabular and FASTA export utilities for MagCluster results."""

from pathlib import Path
from typing import Dict, List
from ..core.cluster import MGCCluster
from ..core.scoring import ClusterScore
from ..core.search import MagnetosomeHit

CORE_FAMILIES = ["MamA", "MamB", "MamE", "MamI", "MamK", "MamL", "MamM", "MamO", "MamP", "MamQ"]

def write_cluster_tables(
    clusters: List[MGCCluster],
    scores: Dict[str, ClusterScore],
    hits: Dict[str, MagnetosomeHit],
    out_dir: Path
):
    """Write standardized TSV reports for clusters and genes."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Clusters table
    clusters_tsv = out_dir / "mgc_clusters.tsv"
    with open(clusters_tsv, "w", encoding="utf-8") as f:
        f.write("cluster_id\tgenome\tcontig\tstart\tend\tspan_bp\tcore_span_bp\tmineral_type\tn_families\tmagnetosome_families\tflags\n")
        for c in clusters:
            s = scores[c.cluster_id]
            f.write(
                f"{c.cluster_id}\t{c.genome_name}\t{c.contig_id}\t{c.start_bp}\t{c.end_bp}\t{c.span_bp}\t"
                f"{c.core_span_bp}\t{s.mineral_type}\t{s.n_mag_families}\t{','.join(s.magnetosome_families_found)}\t{','.join(s.flags)}\n"
            )

    # 2. Genes table
    genes_tsv = out_dir / "mgc_genes.tsv"
    with open(genes_tsv, "w", encoding="utf-8") as f:
        f.write("cluster_id\tgenome\tgene_id\tcontig\tstart\tend\tstrand\tprofile\tparent_family\trole\tscore\tevalue\n")
        for c in clusters:
            for g in c.all_genes:
                if g.gene_id in hits:
                    h = hits[g.gene_id]
                    f.write(
                        f"{c.cluster_id}\t{c.genome_name}\t{g.gene_id}\t{c.contig_id}\t{g.start}\t{g.end}\t"
                        f"{g.strand}\t{h.profile}\t{h.parent_family}\t{h.role}\t{h.score:.1f}\t{h.evalue:.2e}\n"
                    )

    # 3. Presence/Absence Matrix
    matrix_tsv = out_dir / "presence_absence.tsv"
    with open(matrix_tsv, "w", encoding="utf-8") as f:
        f.write("cluster_id\tgenome\t" + "\t".join(CORE_FAMILIES) + "\n")
        for c in clusters:
            s = scores[c.cluster_id]
            vals = [("1" if fam in s.magnetosome_families_found else "0") for fam in CORE_FAMILIES]
            f.write(f"{c.cluster_id}\t{c.genome_name}\t" + "\t".join(vals) + "\n")

    # 4. Summary table (one line per genome)
    summary_tsv = out_dir / "summary.tsv"
    genomes = sorted(list({c.genome_name for c in clusters}))
    with open(summary_tsv, "w", encoding="utf-8") as f:
        f.write("genome\tn_clusters\tmineral_types\n")
        for g in genomes:
            g_clusters = [c for c in clusters if c.genome_name == g]
            g_scores = [scores[c.cluster_id] for c in g_clusters]
            minerals = ",".join(sorted(set(s.mineral_type for s in g_scores)))
            f.write(f"{g}\t{len(g_clusters)}\t{minerals}\n")

def export_combined_magnetosome_fasta(
    clusters: List[MGCCluster],
    hits: Dict[str, MagnetosomeHit],
    out_dir: Path,
    filename: str = "magnetosome_proteins.faa"
):
    """Export all identified magnetosome proteins into a single combined FASTA file."""
    faa_path = out_dir / filename
    seen_genes = set()

    with open(faa_path, "w", encoding="utf-8") as f:
        for c in clusters:
            for g in c.all_genes:
                if g.gene_id in hits and g.gene_id not in seen_genes:
                    seen_genes.add(g.gene_id)
                    h = hits[g.gene_id]
                    hdr = (
                        f">{g.gene_id} [{c.genome_name}] cluster={c.cluster_id} "
                        f"family={h.parent_family} profile={h.profile} score={h.score:.1f} evalue={h.evalue:.2e}"
                    )
                    f.write(f"{hdr}\n{g.translation}\n")
