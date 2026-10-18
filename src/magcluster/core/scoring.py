"""MGC Mineral Type prediction."""

from dataclasses import dataclass
from typing import List
from .cluster import MGCCluster

@dataclass
class ClusterScore:
    """Evaluation and mineral type classification of an MGC."""
    cluster_id: str
    mineral_type: str                   # 'Magnetite', 'Greigite', 'Putative/Mixed'
    magnetosome_families_found: List[str]
    n_mag_families: int
    n_total_genes: int
    flags: List[str]

class ClusterScorer:
    """Annotates MGCs with mineral type prediction."""

    @staticmethod
    def score_cluster(cluster: MGCCluster) -> ClusterScore:
        mag_fams = sorted(list(cluster.distinct_magnetosome_families))
        
        flags = []
        if cluster.on_edge:
            flags.append("Contig_Edge")

        # Mineral type inference
        all_hits = cluster.core_hits + cluster.accessory_hits
        all_profiles = [h.profile.lower() for h in all_hits]
        greigite_hits = sum(1 for p in all_profiles if "greigite" in p)

        if greigite_hits >= 2:
            mineral_type = "Greigite"
        elif greigite_hits == 1:
            mineral_type = "Putative/Mixed"
        else:
            mineral_type = "Magnetite"

        return ClusterScore(
            cluster_id=cluster.cluster_id,
            mineral_type=mineral_type,
            magnetosome_families_found=mag_fams,
            n_mag_families=len(mag_fams),
            n_total_genes=len(cluster.all_genes),
            flags=flags
        )
