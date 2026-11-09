"""Tabular and FASTA export utilities for MagCluster results."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional
from ..core.cluster import MGCCluster
from ..core.predict import PredictedGene
from ..core.scoring import ClusterScore
from ..core.search import MagnetosomeHit


@dataclass
class GenomeSummaryRow:
    """One-row per-genome summary written to ``mgc_summary.tsv``."""
    genome: str
    input_file: str
    status: str
    n_genes: int = 0
    n_hits: int = 0
    n_clustered_hits: int = 0
    cluster_id: str = "-"
    n_subclusters: int = 0
    coordinates: str = "-"
    n_distinct_families: int = 0
    mineral_type: str = "-"
    flags: str = "-"
    families: str = "-"
    gbk_file: str = "-"


SUMMARY_COLUMNS = (
    "genome",
    "input_file",
    "status",
    "n_genes",
    "n_hits",
    "n_clustered_hits",
    "cluster_id",
    "n_subclusters",
    "coordinates",
    "n_distinct_families",
    "mineral_type",
    "flags",
    "families",
    "gbk_file",
)


def write_cohort_summary_table(rows: List[GenomeSummaryRow], out_dir: Path) -> Path:
    """Write ``mgc_summary.tsv`` summarizing every input genome in a batch run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_tsv = out_dir / "mgc_summary.tsv"
    with open(summary_tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(SUMMARY_COLUMNS) + "\n")
        for r in rows:
            vals = [
                r.genome,
                r.input_file,
                r.status,
                str(r.n_genes),
                str(r.n_hits),
                str(r.n_clustered_hits),
                r.cluster_id,
                str(r.n_subclusters),
                r.coordinates,
                str(r.n_distinct_families),
                r.mineral_type,
                r.flags,
                r.families,
                r.gbk_file,
            ]
            f.write("\t".join(vals) + "\n")
    return summary_tsv

def _short_cluster_id(cluster_id: str, genome_name: str) -> str:
    """Strip a redundant leading ``{genome_name}_`` prefix from *cluster_id* for per-genome files."""
    prefix = f"{genome_name}_"
    if genome_name and cluster_id.startswith(prefix) and len(cluster_id) > len(prefix):
        return cluster_id[len(prefix):]
    return cluster_id


def _short_gene_id(gene_id: str, genome_name: str) -> str:
    """Strip a redundant leading ``{genome_name}_`` prefix from *gene_id* for per-genome files."""
    prefix = f"{genome_name}_"
    if genome_name and gene_id.startswith(prefix) and len(gene_id) > len(prefix):
        rest = gene_id[len(prefix):]
        return f"gene_{rest}" if rest.isdigit() else rest
    return gene_id


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
        f.write("cluster_id\tgene_id\tcontig\tstart\tend\tstrand\tprofile\tparent_family\trole\tscore\tevalue\thmm_cov\tseq_cov\n")
        for c in clusters:
            if c.genome_name != genome_name:
                continue
            short_cid = _short_cluster_id(c.cluster_id, genome_name)
            for g in c.all_genes:
                if g.gene_id in hits:
                    written_genes.add(g.gene_id)
                    h = hits[g.gene_id]
                    short_gid = _short_gene_id(g.gene_id, genome_name)
                    strand_str = "+" if g.strand > 0 else ("-" if g.strand < 0 else "-")
                    f.write(
                        f"{short_cid}\t{short_gid}\t{g.contig_id}\t{g.start}\t{g.end}\t"
                        f"{strand_str}\t{h.profile}\t{h.parent_family}\t{h.role}\t{h.score:.1f}\t{h.evalue:.2e}\t"
                        f"{h.hmm_cov:.2f}\t{h.seq_cov:.2f}\n"
                    )

        # Also write unclustered magnetosome hits for this genome if any
        if genes_lookup:
            for gene_id, h in hits.items():
                if gene_id not in written_genes:
                    g = genes_lookup.get(gene_id)
                    if g and getattr(g, "genome_name", "genome") == genome_name:
                        short_gid = _short_gene_id(g.gene_id, genome_name)
                        start_str = str(g.start) if g.start > 0 else "-"
                        end_str = str(g.end) if g.end > 0 else "-"
                        strand_str = "+" if g.strand > 0 else ("-" if g.strand < 0 else "-")
                        f.write(
                            f"unclustered\t{short_gid}\t{g.contig_id}\t{start_str}\t{end_str}\t"
                            f"{strand_str}\t{h.profile}\t{h.parent_family}\t{h.role}\t{h.score:.1f}\t{h.evalue:.2e}\t"
                            f"{h.hmm_cov:.2f}\t{h.seq_cov:.2f}\n"
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
            short_cid = _short_cluster_id(c.cluster_id, genome_name)
            for g in c.all_genes:
                if g.gene_id in hits and g.gene_id not in seen_genes:
                    seen_genes.add(g.gene_id)
                    h = hits[g.gene_id]
                    short_gid = _short_gene_id(g.gene_id, genome_name)
                    hdr = (
                        f">{short_gid} cluster={short_cid} "
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
                        short_gid = _short_gene_id(g.gene_id, genome_name)
                        hdr = (
                            f">{short_gid} cluster=unclustered "
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
