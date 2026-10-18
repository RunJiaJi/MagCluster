"""MGC Completeness, Confidence, Mineral Type, and Lineage classification."""

from collections import Counter
from dataclasses import dataclass
from typing import Dict, List, Set
from .cluster import MGCCluster

TOTAL_CORE_FAMILIES = 10  # MamA, MamB, MamE, MamI, MamK, MamL, MamM, MamO, MamP, MamQ

@dataclass
class ClusterScore:
    """Detailed score and biological classification of an MGC."""
    cluster_id: str
    completeness_score: float           # 0.0 to 1.0 (ratio of detected core families)
    confidence: str                     # 'High', 'Medium', 'Low'
    mineral_type: str                   # 'Magnetite', 'Greigite', 'Mixed/Putative'
    lineage_signature: str              # 'Alphaproteobacteria', 'Desulfobacterota', 'Nitrospirota', 'Environmental_Novel', 'Novel/Mixed'
    core_genes_found: List[str]
    accessory_genes_found: List[str]
    n_core_genes: int
    n_total_genes: int
    flags: List[str]

class ClusterScorer:
    """Evaluates and annotates MGCs with biological metrics."""

    @staticmethod
    def score_cluster(cluster: MGCCluster) -> ClusterScore:
        cores_found = sorted(list(cluster.distinct_core_families))
        acc_found = sorted(list({h.parent_family for h in cluster.accessory_hits}))
        
        n_cores = len(cores_found)
        completeness = min(1.0, n_cores / float(TOTAL_CORE_FAMILIES))
        
        flags = []
        if cluster.on_edge:
            flags.append("Contig_Edge")

        # Confidence assessment
        if n_cores >= 7:
            confidence = "High"
        elif n_cores >= 5:
            confidence = "High" if not cluster.on_edge else "Medium"
        elif n_cores >= 3:
            confidence = "Medium"
        else:
            confidence = "Low"

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

        # Lineage signature: calculate weighted vote across identified genes
        lineage_counts = Counter()
        for h in all_hits:
            lin = h.lineage
            if lin in ("Alpha_classical", "Alpha_Mms"):
                lineage_counts["Alphaproteobacteria"] += 1
            elif lin in ("Desulfo_Mad", "Desulfo_greigite"):
                lineage_counts["Desulfobacterota"] += 1
            elif lin == "Nitrospirota_Man":
                lineage_counts["Nitrospirota"] += 1
            elif lin == "Environmental_Mae":
                lineage_counts["Environmental_Novel"] += 1
            elif lin in ("Novel_Mca", "Novel_Maq"):
                lineage_counts["Novel_Lineage"] += 1

        if lineage_counts:
            # Pick dominant lineage
            dominant_lineage, count = lineage_counts.most_common(1)[0]
            # If close tie between multiple lineages
            total_votes = sum(lineage_counts.values())
            if count / total_votes >= 0.5:
                lineage_sig = dominant_lineage
            else:
                top_two = [l for l, _ in lineage_counts.most_common(2)]
                lineage_sig = f"{top_two[0]}/{top_two[1]}"
        else:
            lineage_sig = "Unknown"

        return ClusterScore(
            cluster_id=cluster.cluster_id,
            completeness_score=round(completeness, 3),
            confidence=confidence,
            mineral_type=mineral_type,
            lineage_signature=lineage_sig,
            core_genes_found=cores_found,
            accessory_genes_found=acc_found,
            n_core_genes=n_cores,
            n_total_genes=len(cluster.all_genes),
            flags=flags
        )
