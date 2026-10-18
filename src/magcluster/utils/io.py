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
    """Write standardized TSV reports for identified magnetosome genes."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # Genes table: mgc_genes.tsv
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
