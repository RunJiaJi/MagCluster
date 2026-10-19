"""Synteny-aware Magnetosome Gene Cluster (MGC) rule engine and clustering."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set
from .predict import PredictedGene
from .search import MagnetosomeHit

@dataclass
class MGCSubcluster:
    """Represents an individual syntenic subcluster (slice) within an MGC."""
    subcluster_id: str
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

    @property
    def span_bp(self) -> int:
        return self.end_bp - self.start_bp + 1

    @property
    def core_span_bp(self) -> int:
        return self.cluster_core_end - self.cluster_core_start + 1

    @property
    def distinct_magnetosome_families(self) -> Set[str]:
        return {h.parent_family for h in (self.core_hits + self.accessory_hits)}

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
    subclusters: List[MGCSubcluster] = field(default_factory=list)
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
        """Group magnetosome hits in a genome into consolidated MGC with discrete subclusters."""
        if not genes or not hits:
            return []

        # Partition genes by contig in genomic order
        genes_by_contig: Dict[str, List[PredictedGene]] = {}
        for g in genes:
            genes_by_contig.setdefault(g.contig_id, []).append(g)

        all_subclusters: List[MGCSubcluster] = []

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

            # Filter groups by minimum distinct magnetosome families
            qualifying_slices = []
            for group in raw_groups:
                first_hit_idx = group[0]
                last_hit_idx = group[-1]

                g_hits = [hits[ctg_genes[idx].gene_id] for idx in group]
                distinct_mag_fams = {h.parent_family for h in g_hits if h.role != "iron_context"}
                if len(distinct_mag_fams) < self.min_distinct_families:
                    continue

                flank_start_idx = max(0, first_hit_idx - self.flank_genes)
                flank_end_idx = min(total_ctg_genes - 1, last_hit_idx + self.flank_genes)

                qualifying_slices.append({
                    "hit_indices": group,
                    "first_hit_idx": first_hit_idx,
                    "last_hit_idx": last_hit_idx,
                    "flank_start_idx": flank_start_idx,
                    "flank_end_idx": flank_end_idx,
                })

            if not qualifying_slices:
                continue

            # Merge contiguous or overlapping flank slices on the same contig
            merged_slices = []
            for q in qualifying_slices:
                if not merged_slices:
                    merged_slices.append(q)
                else:
                    prev = merged_slices[-1]
                    if q["flank_start_idx"] <= prev["flank_end_idx"] + 1:
                        prev["flank_end_idx"] = max(prev["flank_end_idx"], q["flank_end_idx"])
                        prev["last_hit_idx"] = max(prev["last_hit_idx"], q["last_hit_idx"])
                        prev["first_hit_idx"] = min(prev["first_hit_idx"], q["first_hit_idx"])
                        prev["hit_indices"].extend(q["hit_indices"])
                    else:
                        merged_slices.append(q)

            # Build MGCSubcluster objects for each slice
            for s_idx, s in enumerate(merged_slices, 1):
                region_genes = ctg_genes[s["flank_start_idx"] : s["flank_end_idx"] + 1]

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

                on_edge = (s["first_hit_idx"] <= self.flank_genes or s["last_hit_idx"] >= total_ctg_genes - 1 - self.flank_genes)
                start_bp = region_genes[0].start
                end_bp = region_genes[-1].end
                core_start_bp = ctg_genes[s["first_hit_idx"]].start
                core_end_bp = ctg_genes[s["last_hit_idx"]].end

                sub = MGCSubcluster(
                    subcluster_id=f"{genome_name}_loc{s_idx}",
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
                all_subclusters.append(sub)

        if not all_subclusters:
            return []

        # Merge all qualifying subclusters in this genome into exactly ONE consolidated MGC
        cluster_id = f"{genome_name}_MGC"

        consolidated_genes = [g for sub in all_subclusters for g in sub.all_genes]
        consolidated_core = [h for sub in all_subclusters for h in sub.core_hits]
        consolidated_acc = [h for sub in all_subclusters for h in sub.accessory_hits]
        consolidated_iron = [h for sub in all_subclusters for h in sub.iron_hits]
        on_edge = any(sub.on_edge for sub in all_subclusters)

        start_bp = min(sub.start_bp for sub in all_subclusters)
        end_bp = max(sub.end_bp for sub in all_subclusters)
        all_mag_hits = consolidated_core + consolidated_acc
        core_start_bp = min(h.start for h in all_mag_hits) if all_mag_hits else start_bp
        core_end_bp = max(h.end for h in all_mag_hits) if all_mag_hits else end_bp
        contigs_used = list(dict.fromkeys(sub.contig_id for sub in all_subclusters))
        contig_label = ",".join(contigs_used)

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
            subclusters=all_subclusters,
            on_edge=on_edge
        )
        return [final_mgc]
