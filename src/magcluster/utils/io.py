"""Tabular, FASTA, and JSON export utilities for MagCluster results."""

import json
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
        f.write("cluster_id\tgenome\tcontig\tstart\tend\tspan_bp\tcore_span_bp\tcompleteness\tconfidence\tmineral_type\tlineage\tn_core\tcore_genes\taccessory_genes\tflags\n")
        for c in clusters:
            s = scores[c.cluster_id]
            f.write(
                f"{c.cluster_id}\t{c.genome_name}\t{c.contig_id}\t{c.start_bp}\t{c.end_bp}\t{c.span_bp}\t"
                f"{c.core_span_bp}\t{s.completeness_score}\t{s.confidence}\t{s.mineral_type}\t{s.lineage_signature}\t"
                f"{s.n_core_genes}\t{','.join(s.core_genes_found)}\t{','.join(s.accessory_genes_found)}\t{','.join(s.flags)}\n"
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
            vals = [("1" if fam in s.core_genes_found else "0") for fam in CORE_FAMILIES]
            f.write(f"{c.cluster_id}\t{c.genome_name}\t" + "\t".join(vals) + "\n")

    # 4. Summary table (one line per genome)
    summary_tsv = out_dir / "summary.tsv"
    genomes = sorted(list({c.genome_name for c in clusters}))
    with open(summary_tsv, "w", encoding="utf-8") as f:
        f.write("genome\tn_clusters\thighest_confidence\tmineral_types\tlineages\n")
        for g in genomes:
            g_clusters = [c for c in clusters if c.genome_name == g]
            g_scores = [scores[c.cluster_id] for c in g_clusters]
            confidences = [s.confidence for s in g_scores]
            best_conf = "High" if "High" in confidences else ("Medium" if "Medium" in confidences else "Low")
            minerals = ",".join(sorted(set(s.mineral_type for s in g_scores)))
            lineages = ",".join(sorted(set(s.lineage_signature for s in g_scores)))
            f.write(f"{g}\t{len(g_clusters)}\t{best_conf}\t{minerals}\t{lineages}\n")

def export_family_fastas(
    clusters: List[MGCCluster],
    hits: Dict[str, MagnetosomeHit],
    out_dir: Path
):
    """Export extracted protein sequences split by magnetosome family (ready for MSA & phylogenetics)."""
    faa_dir = out_dir / "proteins"
    faa_dir.mkdir(parents=True, exist_ok=True)

    seqs_by_family: Dict[str, List[str]] = {}

    for c in clusters:
        for g in c.all_genes:
            if g.gene_id in hits:
                h = hits[g.gene_id]
                fam = h.parent_family
                hdr = f">{g.gene_id} [{c.genome_name}] cluster={c.cluster_id} family={fam} score={h.score:.1f}"
                seq_entry = f"{hdr}\n{g.translation}\n"
                seqs_by_family.setdefault(fam, []).append(seq_entry)

    for fam, entries in seqs_by_family.items():
        fam_file = faa_dir / f"{fam}.faa"
        with open(fam_file, "w", encoding="utf-8") as f:
            f.writelines(entries)
