"""MGC Mineral Type prediction and quality flag annotation."""

from dataclasses import dataclass, field
from typing import Dict, List, Set
from .cluster import MGCCluster
from .search import MagnetosomeHit

#: Lineage tag identifying greigite-associated HMM profiles in rules.yaml.
GREIGITE_LINEAGE = "Desulfo_greigite"
#: Core families shared by both magnetite and greigite MGCs (no _greigite HMM profile).
SHARED_CORE_FAMILIES = frozenset({"MamK"})
#: Magnetite-specific core families absent from pure greigite MGCs.
MAGNETITE_SPECIFIC_CORE_FAMILIES = frozenset({"MamL", "MamM"})


def _is_greigite_hit(hit: MagnetosomeHit) -> bool:
    """True when a hit matches a greigite-associated HMM profile."""
    return hit.lineage == GREIGITE_LINEAGE or "greigite" in hit.profile.lower()


def _is_complementary_split_only(g_profs: Set[str], r_profs: Set[str]) -> bool:
    """True when greigite and regular hits of a family are strictly complementary Nter/Cter halves.

    Prevents a single split gene (e.g. ``MamE_Nter`` + ``MamE_Cter_greigite``) from being
    misinterpreted as two full dual-mineral copies of the same family.
    """
    if not g_profs or not r_profs:
        return False
    g_nter = all("_nter" in p.lower() for p in g_profs)
    g_cter = all("_cter" in p.lower() for p in g_profs)
    r_nter = all("_nter" in p.lower() for p in r_profs)
    r_cter = all("_cter" in p.lower() for p in r_profs)
    return (g_nter and r_cter) or (g_cter and r_nter)


def _short_subcluster_label(sub_id: str, genome_name: str, fallback_idx: int) -> str:
    """Return concise subcluster label such as 'loc1' for multi-subcluster flags."""
    prefix = f"{genome_name}_"
    if genome_name and sub_id.startswith(prefix) and len(sub_id) > len(prefix):
        return sub_id[len(prefix):]
    return sub_id or f"loc{fallback_idx}"


@dataclass
class ClusterScore:
    """Evaluation and mineral type classification of an MGC."""
    cluster_id: str
    mineral_type: str                   # 'Magnetite', 'Greigite', 'Both', 'Putative_Greigite', 'Unknown'
    magnetosome_families_found: List[str]
    n_mag_families: int
    n_total_genes: int
    flags: List[str]
    greigite_families: List[str] = field(default_factory=list)
    magnetite_core_families: List[str] = field(default_factory=list)
    dual_copy_families: List[str] = field(default_factory=list)


class ClusterScorer:
    """Annotates MGCs with mineral type prediction and structural quality flags."""

    @staticmethod
    def score_cluster(cluster: MGCCluster) -> ClusterScore:
        mag_fams = sorted(list(cluster.distinct_magnetosome_families))

        flags: List[str] = []
        if cluster.subclusters:
            if len(cluster.subclusters) == 1:
                if cluster.subclusters[0].on_edge:
                    flags.append("Contig_Edge")
            else:
                edge_labels = [
                    _short_subcluster_label(sub.subcluster_id, cluster.genome_name, idx)
                    for idx, sub in enumerate(cluster.subclusters, start=1)
                    if sub.on_edge
                ]
                if edge_labels:
                    flags.append(f"Contig_Edge({';'.join(edge_labels)})")
        elif cluster.on_edge:
            flags.append("Contig_Edge")

        if cluster.is_multi_contig:
            flags.append("Multi_Contig")

        if cluster.wraps_origin:
            flags.append("Spans_Origin")

        # Mineral type inference using deduplicated distinct families:
        #   - G: distinct greigite-associated families (Desulfo_greigite, incl. MamT_greigite)
        #   - M: distinct non-greigite core families excluding shared MamK and complementary split halves
        #   - M_all: distinct non-greigite core families including MamK (10 families)
        #   - dual_copy_fams: core families present in both regular and _greigite copies
        #   - magnetite_specific: MamL / MamM (absent from pure greigite MGCs)
        all_hits = cluster.core_hits + cluster.accessory_hits
        greigite_profs_by_fam: Dict[str, Set[str]] = {}
        regular_core_profs_by_fam: Dict[str, Set[str]] = {}
        for h in all_hits:
            if _is_greigite_hit(h):
                greigite_profs_by_fam.setdefault(h.parent_family, set()).add(h.profile)
            elif h.role == "core":
                regular_core_profs_by_fam.setdefault(h.parent_family, set()).add(h.profile)

        greigite_fam_set = set(greigite_profs_by_fam.keys())
        non_greigite_all_core_fams = set(regular_core_profs_by_fam.keys())
        shared_core_fams = greigite_fam_set & non_greigite_all_core_fams
        split_only_fams = {
            fam
            for fam in shared_core_fams
            if _is_complementary_split_only(
                greigite_profs_by_fam[fam], regular_core_profs_by_fam[fam]
            )
        }
        # Exclude shared MamK and complementary split-gene halves (e.g. MamE_Nter + MamE_Cter_greigite)
        # so one split gene cannot double-count into both G and M/dual_copy.
        magnetite_core_fam_set = non_greigite_all_core_fams - SHARED_CORE_FAMILIES - split_only_fams
        dual_copy_fam_set = shared_core_fams - split_only_fams
        has_mag_specific_core = bool(magnetite_core_fam_set & MAGNETITE_SPECIFIC_CORE_FAMILIES)

        n_greigite = len(greigite_fam_set)
        n_mag_core = len(magnetite_core_fam_set)
        n_all_non_greigite_core = len(non_greigite_all_core_fams)

        if n_greigite >= 3:
            if (bool(dual_copy_fam_set) and n_mag_core >= 2) or (
                has_mag_specific_core and n_mag_core >= 3
            ):
                mineral_type = "Both"
            else:
                mineral_type = "Greigite"
        elif n_greigite >= 1:
            if n_mag_core >= 3:
                mineral_type = "Magnetite"
            else:
                mineral_type = "Putative_Greigite"
        else:
            if n_all_non_greigite_core >= 3:
                mineral_type = "Magnetite"
            else:
                mineral_type = "Unknown"

        return ClusterScore(
            cluster_id=cluster.cluster_id,
            mineral_type=mineral_type,
            magnetosome_families_found=mag_fams,
            n_mag_families=len(mag_fams),
            n_total_genes=len(cluster.all_genes),
            flags=flags,
            greigite_families=sorted(greigite_fam_set),
            magnetite_core_families=sorted(magnetite_core_fam_set),
            dual_copy_families=sorted(dual_copy_fam_set),
        )
