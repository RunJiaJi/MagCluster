"""Tabular and FASTA export utilities for MagCluster results."""

from pathlib import Path
from typing import Dict, List, Optional
from ..core.cluster import MGCCluster
from ..core.predict import PredictedGene
from ..core.scoring import ClusterScore
from ..core.search import MagnetosomeHit

CORE_FAMILIES = ["MamA", "MamB", "MamE", "MamI", "MamK", "MamL", "MamM", "MamO", "MamP", "MamQ"]

def write_genome_cluster_tables(
    genome_name: str,
    clusters: List[MGCCluster],
    hits: Dict[str, MagnetosomeHit],
    out_dir: Path,
    genes_lookup: Optional[Dict[str, PredictedGene]] = None
) -> Path:
    """Write standardized TSV report for identified magnetosome genes in a specific genome."""
    out_dir.mkdir(parents=True, exist_ok=True)

    genes_tsv = out_dir / f"{genome_name}_mgc_genes.tsv"
    written_genes = set()

    with open(genes_tsv, "w", encoding="utf-8") as f:
        f.write("cluster_id\tgenome\tgene_id\tcontig\tstart\tend\tstrand\tprofile\tparent_family\trole\tscore\tevalue\n")
        for c in clusters:
            if c.genome_name != genome_name:
                continue
            for g in c.all_genes:
                if g.gene_id in hits:
                    written_genes.add(g.gene_id)
                    h = hits[g.gene_id]
                    strand_str = "+" if g.strand > 0 else ("-" if g.strand < 0 else "-")
                    f.write(
                        f"{c.cluster_id}\t{c.genome_name}\t{g.gene_id}\t{g.contig_id}\t{g.start}\t{g.end}\t"
                        f"{strand_str}\t{h.profile}\t{h.parent_family}\t{h.role}\t{h.score:.1f}\t{h.evalue:.2e}\n"
                    )

        # Also write unclustered magnetosome hits for this genome if any
        if genes_lookup:
            for gene_id, h in hits.items():
                if gene_id not in written_genes:
                    g = genes_lookup.get(gene_id)
                    if g and getattr(g, "genome_name", "genome") == genome_name:
                        start_str = str(g.start) if g.start > 0 else "-"
                        end_str = str(g.end) if g.end > 0 else "-"
                        strand_str = "+" if g.strand > 0 else ("-" if g.strand < 0 else "-")
                        f.write(
                            f"unclustered\t{genome_name}\t{g.gene_id}\t{g.contig_id}\t{start_str}\t{end_str}\t"
                            f"{strand_str}\t{h.profile}\t{h.parent_family}\t{h.role}\t{h.score:.1f}\t{h.evalue:.2e}\n"
                        )
    return genes_tsv

def export_genome_magnetosome_fasta(
    genome_name: str,
    clusters: List[MGCCluster],
    hits: Dict[str, MagnetosomeHit],
    out_dir: Path,
    genes_lookup: Optional[Dict[str, PredictedGene]] = None
) -> Path:
    """Export identified magnetosome proteins for a specific genome into FASTA format."""
    out_dir.mkdir(parents=True, exist_ok=True)
    faa_path = out_dir / f"{genome_name}_magnetosome_proteins.faa"
    seen_genes = set()

    with open(faa_path, "w", encoding="utf-8") as f:
        for c in clusters:
            if c.genome_name != genome_name:
                continue
            for g in c.all_genes:
                if g.gene_id in hits and g.gene_id not in seen_genes:
                    seen_genes.add(g.gene_id)
                    h = hits[g.gene_id]
                    hdr = (
                        f">{g.gene_id} [{c.genome_name}] cluster={c.cluster_id} "
                        f"family={h.parent_family} profile={h.profile} score={h.score:.1f} evalue={h.evalue:.2e}"
                    )
                    f.write(f"{hdr}\n{g.translation}\n")

        # Also export unclustered magnetosome hits for this genome if any
        if genes_lookup:
            for gene_id, h in hits.items():
                if gene_id not in seen_genes:
                    g = genes_lookup.get(gene_id)
                    if g and getattr(g, "genome_name", "genome") == genome_name:
                        seen_genes.add(gene_id)
                        hdr = (
                            f">{g.gene_id} [{genome_name}] cluster=unclustered "
                            f"family={h.parent_family} profile={h.profile} score={h.score:.1f} evalue={h.evalue:.2e}"
                        )
                        f.write(f"{hdr}\n{g.translation}\n")
    return faa_path

def write_all_genome_tables(
    genome_names: List[str],
    clusters: List[MGCCluster],
    scores: Dict[str, ClusterScore],
    hits: Dict[str, MagnetosomeHit],
    out_dir: Path,
    genes_lookup: Optional[Dict[str, PredictedGene]] = None
) -> List[Path]:
    """Write per-genome standardized TSV reports for identified magnetosome genes."""
    out_paths = []
    for gname in genome_names:
        p = write_genome_cluster_tables(gname, clusters, hits, out_dir, genes_lookup=genes_lookup)
        out_paths.append(p)
    return out_paths

def export_all_genome_fastas(
    genome_names: List[str],
    clusters: List[MGCCluster],
    hits: Dict[str, MagnetosomeHit],
    out_dir: Path,
    genes_lookup: Optional[Dict[str, PredictedGene]] = None
) -> List[Path]:
    """Export per-genome identified magnetosome proteins into separate FASTA files."""
    out_paths = []
    for gname in genome_names:
        p = export_genome_magnetosome_fasta(gname, clusters, hits, out_dir, genes_lookup=genes_lookup)
        out_paths.append(p)
    return out_paths

# Backward compatibility wrappers
def write_cluster_tables(
    clusters: List[MGCCluster],
    scores: Dict[str, ClusterScore],
    hits: Dict[str, MagnetosomeHit],
    out_dir: Path,
    genes_lookup: Optional[Dict[str, PredictedGene]] = None
):
    genomes = list(dict.fromkeys([c.genome_name for c in clusters]))
    if not genomes and genes_lookup:
        genomes = list(dict.fromkeys([getattr(g, "genome_name", "genome") for g in genes_lookup.values() if g.gene_id in hits]))
    write_all_genome_tables(genomes, clusters, scores, hits, out_dir, genes_lookup=genes_lookup)

def export_combined_magnetosome_fasta(
    clusters: List[MGCCluster],
    hits: Dict[str, MagnetosomeHit],
    out_dir: Path,
    filename: str = "magnetosome_proteins.faa",
    genes_lookup: Optional[Dict[str, PredictedGene]] = None
):
    genomes = list(dict.fromkeys([c.genome_name for c in clusters]))
    if not genomes and genes_lookup:
        genomes = list(dict.fromkeys([getattr(g, "genome_name", "genome") for g in genes_lookup.values() if g.gene_id in hits]))
    export_all_genome_fastas(genomes, clusters, hits, out_dir, genes_lookup=genes_lookup)
