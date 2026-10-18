"""Synteny-aware Magnetosome Gene Cluster (MGC) rule engine and clustering."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set
from .predict import PredictedGene
from .search import MagnetosomeHit

@dataclass
class MGCCluster:
    """Represents an identified Magnetosome Gene Cluster."""
    cluster_id: str
    genome_name: str
    contig_id: str
    start_bp: int                       # 1-based start (including flanks)
    end_bp: int                         # 1-based end (including flanks)
    cluster_core_start: int             # start of first magnetosome gene
    cluster_core_end: int               # end of last magnetosome gene
    core_hits: List[MagnetosomeHit] = field(default_factory=list)
    accessory_hits: List[MagnetosomeHit] = field(default_factory=list)
    iron_hits: List[MagnetosomeHit] = field(default_factory=list)
    all_genes: List[PredictedGene] = field(default_factory=list)
    on_edge: bool = False
    is_circular: bool = False

    @property
    def distinct_magnetosome_families(self) -> Set[str]:
        """Set of distinct true magnetosome families (excluding iron_context/Feo)."""
        return {h.parent_family for h in (self.core_hits + self.accessory_hits)}

    @property
    def distinct_all_families(self) -> Set[str]:
        return {h.parent_family for h in (self.core_hits + self.accessory_hits + self.iron_hits)}

    @property
    def span_bp(self) -> int:
        return self.end_bp - self.start_bp + 1

    @property
    def core_span_bp(self) -> int:
        return self.cluster_core_end - self.cluster_core_start + 1

class ClusterEngine:
    """Detects MGCs based on synteny, gene distances, and magnetosome family threshold."""

    def __init__(
        self,
        max_intergenic_genes: int = 5,
        max_distance_bp: int = 25_000,
        min_distinct_families: int = 3,
        flank_genes: int = 3
    ):
        self.max_intergenic_genes = max_intergenic_genes
        self.max_distance_bp = max_distance_bp
        self.min_distinct_families = min_distinct_families
        self.flank_genes = flank_genes

    def cluster_genome(
        self,
        genome_name: str,
        genes: List[PredictedGene],
        hits: Dict[str, MagnetosomeHit]
    ) -> List[MGCCluster]:
        """Group magnetosome hits in a genome into exactly one consolidated MGC per genome."""
        if not genes or not hits:
            return []

        # Partition genes by contig in genomic order
        genes_by_contig: Dict[str, List[PredictedGene]] = {}
        for g in genes:
            genes_by_contig.setdefault(g.contig_id, []).append(g)

        raw_clusters: List[MGCCluster] = []

        for ctg_id, ctg_genes in genes_by_contig.items():
            ctg_genes.sort(key=lambda x: x.start)
            total_ctg_genes = len(ctg_genes)

            # Find indices of genes with true magnetosome hits (excluding lone Feo)
            hit_indices = []
            for idx, g in enumerate(ctg_genes):
                if g.gene_id in hits:
                    hit = hits[g.gene_id]
                    if hit.role != "iron_context":
                        hit_indices.append(idx)

            if not hit_indices:
                continue

            # Group hit indices by synteny (max_intergenic_genes & max_distance_bp)
            raw_groups: List[List[int]] = []
            curr_group = [hit_indices[0]]

            for prev_idx, next_idx in zip(hit_indices[:-1], hit_indices[1:]):
                prev_g = ctg_genes[prev_idx]
                next_g = ctg_genes[next_idx]
                gene_gap = next_idx - prev_idx - 1
                bp_gap = next_g.start - prev_g.end

                if gene_gap <= self.max_intergenic_genes and bp_gap <= self.max_distance_bp:
                    curr_group.append(next_idx)
                else:
                    raw_groups.append(curr_group)
                    curr_group = [next_idx]
            raw_groups.append(curr_group)

            for group in raw_groups:
                first_hit_idx = group[0]
                last_hit_idx = group[-1]

                flank_start_idx = max(0, first_hit_idx - self.flank_genes)
                flank_end_idx = min(total_ctg_genes - 1, last_hit_idx + self.flank_genes)

                region_genes = ctg_genes[flank_start_idx : flank_end_idx + 1]

                core_hits = []
                acc_hits = []
                iron_hits = []

                for g in region_genes:
                    if g.gene_id in hits:
                        h = hits[g.gene_id]
                        if h.role == "core":
                            core_hits.append(h)
                        elif h.role == "accessory":
                            acc_hits.append(h)
                        elif h.role == "iron_context":
                            iron_hits.append(h)

                distinct_mag_fams = {h.parent_family for h in (core_hits + acc_hits)}
                if len(distinct_mag_fams) < self.min_distinct_families:
                    continue

                on_edge = (first_hit_idx <= self.flank_genes or last_hit_idx >= total_ctg_genes - 1 - self.flank_genes)
                start_bp = region_genes[0].start
                end_bp = region_genes[-1].end
                core_start_bp = ctg_genes[first_hit_idx].start
                core_end_bp = ctg_genes[last_hit_idx].end

                mgc = MGCCluster(
                    cluster_id="tmp",
                    genome_name=genome_name,
                    contig_id=ctg_id,
                    start_bp=start_bp,
                    end_bp=end_bp,
                    cluster_core_start=core_start_bp,
                    cluster_core_end=core_end_bp,
                    core_hits=core_hits,
                    accessory_hits=acc_hits,
                    iron_hits=iron_hits,
                    all_genes=region_genes,
                    on_edge=on_edge
                )
                raw_clusters.append(mgc)

        if not raw_clusters:
            return []

        # Merge all qualifying sub-clusters in this genome into exactly ONE consolidated MGC
        cluster_id = f"{genome_name}_MGC"

        contigs_used = []
        consolidated_genes: List[PredictedGene] = []
        consolidated_core: List[MagnetosomeHit] = []
        consolidated_acc: List[MagnetosomeHit] = []
        consolidated_iron: List[MagnetosomeHit] = []
        on_edge = any(c.on_edge for c in raw_clusters)

        for ctg_id, ctg_genes in genes_by_contig.items():
            ctg_sub_clusters = [c for c in raw_clusters if c.contig_id == ctg_id]
            if not ctg_sub_clusters:
                continue

            contigs_used.append(ctg_id)
            min_start = min(c.start_bp for c in ctg_sub_clusters)
            max_end = max(c.end_bp for c in ctg_sub_clusters)

            # Get all genes on this contig within [min_start, max_end]
            span_genes = [g for g in ctg_genes if g.start >= min_start and g.end <= max_end]
            consolidated_genes.extend(span_genes)

            for g in span_genes:
                if g.gene_id in hits:
                    h = hits[g.gene_id]
                    if h.role == "core":
                        consolidated_core.append(h)
                    elif h.role == "accessory":
                        consolidated_acc.append(h)
                    elif h.role == "iron_context":
                        consolidated_iron.append(h)

        start_bp = min(g.start for g in consolidated_genes)
        end_bp = max(g.end for g in consolidated_genes)
        all_mag_hits = consolidated_core + consolidated_acc
        core_start_bp = min(h.start for h in all_mag_hits) if all_mag_hits else start_bp
        core_end_bp = max(h.end for h in all_mag_hits) if all_mag_hits else end_bp
        contig_label = ",".join(contigs_used) if len(contigs_used) > 1 else contigs_used[0]

        final_mgc = MGCCluster(
            cluster_id=cluster_id,
            genome_name=genome_name,
            contig_id=contig_label,
            start_bp=start_bp,
            end_bp=end_bp,
            cluster_core_start=core_start_bp,
            cluster_core_end=core_end_bp,
            core_hits=consolidated_core,
            accessory_hits=consolidated_acc,
            iron_hits=consolidated_iron,
            all_genes=consolidated_genes,
            on_edge=on_edge
        )
        return [final_mgc]
