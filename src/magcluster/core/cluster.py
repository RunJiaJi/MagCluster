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

    @property
    def clustered_hit_ids(self) -> Set[str]:
        """Gene IDs of HMM hits that belong to this subcluster (excluding unclustered flank hits)."""
        return {h.gene_id for h in (self.core_hits + self.accessory_hits + self.iron_hits)}

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
    #: True when at least one circular replicon has boundary subclusters bridging
    #: across the 0 bp origin.
    wraps_origin: bool = False
    #: True when ``wraps_origin`` is True and *every* subcluster lies within one compact
    #: origin-spanning arc (``start_bp > end_bp``). False if additional distant islands exist.
    wraps_origin_compact: bool = False
    #: Length of the replicon (only set when ``wraps_origin`` is True on a single-contig MGC).
    contig_length: int = 0
    #: Per-contig replicon lengths for circular contigs whose subclusters bridge the 0 bp origin.
    origin_contigs: Dict[str, int] = field(default_factory=dict)

    @property
    def distinct_magnetosome_families(self) -> Set[str]:
        """Set of distinct true magnetosome families (excluding iron_context/Feo)."""
        return {h.parent_family for h in (self.core_hits + self.accessory_hits)}

    @property
    def distinct_all_families(self) -> Set[str]:
        return {h.parent_family for h in (self.core_hits + self.accessory_hits + self.iron_hits)}

    @property
    def distinct_contigs(self) -> List[str]:
        """Ordered list of distinct contig IDs spanned by this MGC."""
        if self.subclusters:
            return list(dict.fromkeys(sub.contig_id for sub in self.subclusters))
        return [c for c in self.contig_id.split(",") if c]

    @property
    def is_multi_contig(self) -> bool:
        """True when MGC subclusters span more than one contig."""
        return len(self.distinct_contigs) > 1

    @property
    def edge_subclusters(self) -> List[MGCSubcluster]:
        """Subclusters whose magnetosome gene sits directly at a linear contig boundary."""
        return [sub for sub in self.subclusters if sub.on_edge]

    @property
    def clustered_hit_ids(self) -> Set[str]:
        """Gene IDs of HMM hits belonging to this MGC (excluding unclustered flank hits)."""
        return {h.gene_id for h in (self.core_hits + self.accessory_hits + self.iron_hits)}

    @property
    def span_bp(self) -> int:
        if self.wraps_origin and self.wraps_origin_compact and self.contig_length:
            return self.contig_length - self.start_bp + 1 + self.end_bp
        if self.wraps_origin and self.subclusters:
            if self.is_multi_contig and self.origin_contigs:
                by_ctg: Dict[str, List[MGCSubcluster]] = {}
                for sub in self.subclusters:
                    by_ctg.setdefault(sub.contig_id, []).append(sub)
                total = 0
                for cid, csubs in by_ctg.items():
                    clen = self.origin_contigs.get(cid, 0)
                    if clen and len(csubs) >= 2:
                        ordered = sorted(csubs, key=lambda x: x.start_bp)
                        total += clen - ordered[-1].start_bp + 1 + ordered[0].end_bp
                        total += sum(s.span_bp for s in ordered[1:-1])
                    else:
                        total += sum(s.span_bp for s in csubs)
                return total
            if self.contig_length and len(self.subclusters) >= 3:
                subs = sorted(self.subclusters, key=lambda x: x.start_bp)
                origin_arc = self.contig_length - subs[-1].start_bp + 1 + subs[0].end_bp
                return origin_arc + sum(sub.span_bp for sub in subs[1:-1])
        if self.subclusters and len(self.subclusters) > 1:
            return sum(sub.span_bp for sub in self.subclusters)
        return self.end_bp - self.start_bp + 1

    @property
    def core_span_bp(self) -> int:
        if (
            self.wraps_origin
            and self.wraps_origin_compact
            and self.contig_length
            and self.cluster_core_start > self.cluster_core_end
        ):
            return self.contig_length - self.cluster_core_start + 1 + self.cluster_core_end
        if self.wraps_origin and self.subclusters:
            if self.is_multi_contig and self.origin_contigs:
                by_ctg: Dict[str, List[MGCSubcluster]] = {}
                for sub in self.subclusters:
                    by_ctg.setdefault(sub.contig_id, []).append(sub)
                total = 0
                for cid, csubs in by_ctg.items():
                    clen = self.origin_contigs.get(cid, 0)
                    if clen and len(csubs) >= 2:
                        ordered = sorted(csubs, key=lambda x: x.start_bp)
                        total += clen - ordered[-1].cluster_core_start + 1 + ordered[0].cluster_core_end
                        total += sum(s.core_span_bp for s in ordered[1:-1])
                    else:
                        total += sum(s.core_span_bp for s in csubs)
                return total
            if self.contig_length and len(self.subclusters) >= 3:
                subs = sorted(self.subclusters, key=lambda x: x.start_bp)
                origin_core = self.contig_length - subs[-1].cluster_core_start + 1 + subs[0].cluster_core_end
                return origin_core + sum(sub.core_span_bp for sub in subs[1:-1])
        if self.subclusters and len(self.subclusters) > 1:
            return sum(sub.core_span_bp for sub in self.subclusters)
        return self.cluster_core_end - self.cluster_core_start + 1

    def format_coordinates(self, sep: str = ", ") -> str:
        """Format cluster genomic coordinates, listing discrete subcluster intervals when > 1."""
        if self.wraps_origin and self.subclusters and self.is_multi_contig and self.origin_contigs:
            by_ctg: Dict[str, List[MGCSubcluster]] = {}
            for sub in self.subclusters:
                by_ctg.setdefault(sub.contig_id, []).append(sub)
            parts: List[str] = []
            for cid, csubs in by_ctg.items():
                if cid in self.origin_contigs and len(csubs) >= 2:
                    ordered = sorted(csubs, key=lambda x: x.start_bp)
                    parts.extend(f"{s.contig_id}:{s.start_bp}-{s.end_bp}" for s in ordered[1:-1])
                    parts.append(f"{ordered[-1].contig_id}:{ordered[-1].start_bp}-{ordered[0].end_bp} (spans origin)")
                else:
                    parts.extend(f"{s.contig_id}:{s.start_bp}-{s.end_bp}" for s in csubs)
            return sep.join(parts)
        if self.wraps_origin and self.subclusters and len(self.subclusters) >= 3 and not self.wraps_origin_compact:
            subs = sorted(self.subclusters, key=lambda x: x.start_bp)
            parts = [f"{sub.contig_id}:{sub.start_bp}-{sub.end_bp}" for sub in subs[1:-1]]
            parts.append(f"{subs[-1].contig_id}:{subs[-1].start_bp}-{subs[0].end_bp} (spans origin)")
            return sep.join(parts)
        if not self.wraps_origin and self.subclusters and len(self.subclusters) > 1:
            return sep.join(f"{sub.contig_id}:{sub.start_bp}-{sub.end_bp}" for sub in self.subclusters)
        origin = " (spans origin)" if self.wraps_origin else ""
        return f"{self.contig_id}:{self.start_bp}-{self.end_bp}{origin}"

def _origin_wrapped_span(
    subclusters: List[MGCSubcluster],
    contig_length: int,
    max_origin_gap: Optional[int] = None,
    force_origin_bridge: bool = False,
):
    """Detect whether *subclusters* on a circular replicon bridge across the 0 bp origin.

    Returns ``(start_bp, end_bp, core_start_bp, core_end_bp, is_compact)`` or ``None``.
    - ``is_compact=True`` when *all* subclusters lie within a single origin-spanning arc
      (every non-largest internal gap is also ``<= max_origin_gap``).
    - ``is_compact=False`` when ``len(subclusters) >= 3`` and the boundary pair bridges
      across the origin while one or more distant independent islands lie elsewhere on the
      chromosome (so distant islands are NOT collapsed into a megabase arc).
    - ``force_origin_bridge=True`` marks replicons whose boundary groups were syntenically
      rescued across the 0 bp origin (preserving ``Spans_Origin`` even on small plasmids
      where the non-MGC interior bp span happens to be shorter than ``origin_gap``).
    """
    if contig_length <= 0 or len(subclusters) < 2:
        return None
    subs = sorted(subclusters, key=lambda x: x.cluster_core_start)
    internal_gaps = [
        subs[i + 1].cluster_core_start - subs[i].cluster_core_end - 1
        for i in range(len(subs) - 1)
    ]
    best_gap = max(internal_gaps)
    best_i = internal_gaps.index(best_gap)
    origin_gap = (
        (contig_length - max(s.cluster_core_end for s in subs))
        + (min(s.cluster_core_start for s in subs) - 1)
    )
    if (not force_origin_bridge and best_gap <= origin_gap) or (
        max_origin_gap is not None and origin_gap > max_origin_gap
    ):
        return None

    other_gaps = [g for idx, g in enumerate(internal_gaps) if idx != best_i]
    is_compact = (
        max_origin_gap is None
        or all(g <= max_origin_gap for g in other_gaps)
    )
    if is_compact:
        return (
            subs[best_i + 1].start_bp,
            subs[best_i].end_bp,
            subs[best_i + 1].cluster_core_start,
            subs[best_i].cluster_core_end,
            True,
        )
    return (
        subs[-1].start_bp,
        subs[0].end_bp,
        subs[-1].cluster_core_start,
        subs[0].cluster_core_end,
        False,
    )

class ClusterEngine:
    """Detects MGCs based on synteny, gene distances, and magnetosome family threshold."""

    def __init__(
        self,
        max_intergenic_genes: int = 5,
        max_distance_bp: int = 25_000,
        min_distinct_families: int = 3,
        flank_genes: int = 3
    ):
        if max_intergenic_genes < 0:
            raise ValueError("max_intergenic_genes must be >= 0")
        if max_distance_bp < 0:
            raise ValueError("max_distance_bp must be >= 0")
        if min_distinct_families < 1:
            raise ValueError("min_distinct_families must be >= 1")
        if flank_genes < 0:
            raise ValueError("flank_genes must be >= 0")
        self.max_intergenic_genes = max_intergenic_genes
        self.max_distance_bp = max_distance_bp
        self.min_distinct_families = min_distinct_families
        self.flank_genes = flank_genes

    def cluster_genome(
        self,
        genome_name: str,
        genes: List[PredictedGene],
        hits: Dict[str, MagnetosomeHit],
        contig_lengths: Optional[Dict[str, int]] = None,
        circular_contigs: Optional[Dict[str, bool]] = None,
    ) -> List[MGCCluster]:
        """Group magnetosome hits in a genome into consolidated MGC with discrete subclusters.

        ``contig_lengths``/``circular_contigs`` are optional and only used to report
        origin-spanning coordinates on complete circular replicons.
        """
        if not genes or not hits:
            return []

        # Partition genes by contig in genomic order (skip coordinate-less genes, e.g. mixed .faa)
        genes_by_contig: Dict[str, List[PredictedGene]] = {}
        for g in genes:
            if g.start <= 0 and g.end <= 0:
                continue
            genes_by_contig.setdefault(g.contig_id, []).append(g)

        all_subclusters: List[MGCSubcluster] = []
        origin_bridged_contigs: Set[str] = set()

        for ctg_id, ctg_genes in genes_by_contig.items():
            ctg_genes.sort(key=lambda x: (x.start, x.end))
            total_ctg_genes = len(ctg_genes)
            is_ctg_circular = bool(circular_contigs and circular_contigs.get(ctg_id))
            ctg_len = (contig_lengths or {}).get(ctg_id, 0)

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
                bp_gap = next_g.start - prev_g.end - 1

                if gene_gap <= self.max_intergenic_genes and bp_gap <= self.max_distance_bp:
                    curr_group.append(next_idx)
                else:
                    raw_groups.append(curr_group)
                    curr_group = [next_idx]
            raw_groups.append(curr_group)

            # On a complete circular replicon, check whether the first and last linear groups
            # are syntenically contiguous across the 0 bp origin (rescuing MGCs split across origin
            # where neither side alone meets min_distinct_families, while their union does).
            origin_bridged_groups = set()
            if is_ctg_circular and ctg_len > 0 and len(raw_groups) >= 2:
                first_grp, last_grp = raw_groups[0], raw_groups[-1]
                origin_gene_gap = (total_ctg_genes - 1 - last_grp[-1]) + first_grp[0]
                last_end = max(ctg_genes[i].end for i in last_grp)
                first_start = min(ctg_genes[i].start for i in first_grp)
                origin_bp_gap = (ctg_len - last_end) + (first_start - 1)
                if origin_gene_gap <= self.max_intergenic_genes and origin_bp_gap <= self.max_distance_bp:
                    combined_fams = {
                        hits[ctg_genes[i].gene_id].parent_family
                        for i in (first_grp + last_grp)
                        if hits[ctg_genes[i].gene_id].role != "iron_context"
                    }
                    if len(combined_fams) >= self.min_distinct_families:
                        origin_bridged_groups = {0, len(raw_groups) - 1}
                        origin_bridged_contigs.add(ctg_id)

            # Filter groups by minimum distinct magnetosome families
            qualifying_slices = []
            for grp_idx, group in enumerate(raw_groups):
                first_hit_idx = group[0]
                last_hit_idx = group[-1]

                g_hits = [hits[ctg_genes[idx].gene_id] for idx in group]
                distinct_mag_fams = {h.parent_family for h in g_hits if h.role != "iron_context"}
                if grp_idx not in origin_bridged_groups and len(distinct_mag_fams) < self.min_distinct_families:
                    continue

                # Extend up to flank_genes on each side, stopping if an inter-gene gap
                # exceeds max_distance_bp (e.g. a large repeat or unannotated scaffold N-gap).
                flank_start_idx = first_hit_idx
                for _ in range(self.flank_genes):
                    if flank_start_idx <= 0:
                        break
                    if ctg_genes[flank_start_idx].start - ctg_genes[flank_start_idx - 1].end - 1 > self.max_distance_bp:
                        break
                    flank_start_idx -= 1

                flank_end_idx = last_hit_idx
                for _ in range(self.flank_genes):
                    if flank_end_idx >= total_ctg_genes - 1:
                        break
                    if ctg_genes[flank_end_idx + 1].start - ctg_genes[flank_end_idx].end - 1 > self.max_distance_bp:
                        break
                    flank_end_idx += 1

                qualifying_slices.append({
                    "hit_indices": list(group),
                    "first_hit_idx": first_hit_idx,
                    "last_hit_idx": last_hit_idx,
                    "flank_start_idx": flank_start_idx,
                    "flank_end_idx": flank_end_idx,
                    "origin_bridged": grp_idx in origin_bridged_groups,
                })

            if not qualifying_slices:
                continue

            # Merge contiguous or overlapping flank slices on the same contig
            # (never bridge slices separated by more than max_distance_bp, and never
            # linearly merge two origin-bridged boundary slices across the non-MGC interior).
            merged_slices = []
            for q in qualifying_slices:
                if not merged_slices:
                    merged_slices.append(q)
                else:
                    prev = merged_slices[-1]
                    prev_max_end = max(
                        ctg_genes[i].end
                        for i in range(prev["flank_start_idx"], prev["flank_end_idx"] + 1)
                    )
                    slice_bp_gap = ctg_genes[q["flank_start_idx"]].start - prev_max_end - 1
                    both_origin_bridged = prev.get("origin_bridged") and q.get("origin_bridged")
                    if (
                        not both_origin_bridged
                        and q["flank_start_idx"] <= prev["flank_end_idx"] + 1
                        and slice_bp_gap <= self.max_distance_bp
                    ):
                        prev["flank_end_idx"] = max(prev["flank_end_idx"], q["flank_end_idx"])
                        prev["last_hit_idx"] = max(prev["last_hit_idx"], q["last_hit_idx"])
                        prev["first_hit_idx"] = min(prev["first_hit_idx"], q["first_hit_idx"])
                        prev["hit_indices"].extend(q["hit_indices"])
                        prev["origin_bridged"] = bool(prev.get("origin_bridged") or q.get("origin_bridged"))
                    else:
                        if both_origin_bridged and q["flank_start_idx"] <= prev["flank_end_idx"]:
                            # Partition interior flanking genes cleanly so origin-bridged
                            # slices remain disjoint (no duplicated genes or span > contig_len).
                            mid = (prev["last_hit_idx"] + q["first_hit_idx"]) // 2
                            prev["flank_end_idx"] = min(prev["flank_end_idx"], mid)
                            q["flank_start_idx"] = max(q["flank_start_idx"], mid + 1)
                        merged_slices.append(q)

            # Build MGCSubcluster objects for each slice (numbered sequentially across the whole genome)
            for s in merged_slices:
                sub_idx = len(all_subclusters) + 1
                region_genes = ctg_genes[s["flank_start_idx"] : s["flank_end_idx"] + 1]
                hit_idx_set = set(s["hit_indices"])
                core_first = s["first_hit_idx"]
                core_last = s["last_hit_idx"]

                core_hits = []
                acc_hits = []
                iron_hits = []

                # Only hits that passed synteny clustering (hit_indices) belong to core/accessory hits;
                # iron_context (Feo) hits inside the core span are included, while distant unclustered
                # hits that merely fall inside --flank remain unclustered context.
                for ctg_idx in range(s["flank_start_idx"], s["flank_end_idx"] + 1):
                    g = ctg_genes[ctg_idx]
                    if g.gene_id not in hits:
                        continue
                    h = hits[g.gene_id]
                    if ctg_idx in hit_idx_set:
                        if h.role == "core":
                            core_hits.append(h)
                        elif h.role == "accessory":
                            acc_hits.append(h)
                    elif h.role == "iron_context" and core_first <= ctg_idx <= core_last:
                        iron_hits.append(h)

                # Contig_Edge is decoupled from --flank and applies to linear contigs only:
                # flagged iff a clustered magnetosome protein sits directly at the first or
                # last gene position on that contig (indicating potential truncation).
                on_edge = (
                    not is_ctg_circular
                    and (core_first == 0 or core_last == total_ctg_genes - 1)
                )
                start_bp = min(g.start for g in region_genes)
                end_bp = max(g.end for g in region_genes)
                core_start_bp = min(ctg_genes[i].start for i in s["hit_indices"])
                core_end_bp = max(ctg_genes[i].end for i in s["hit_indices"])

                sub = MGCSubcluster(
                    subcluster_id=f"{genome_name}_loc{sub_idx}",
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
        # Avoid coarse genome-wide any() aggregation for multi-subcluster MGCs:
        # cluster.on_edge is only set for single-subcluster MGCs, while multi-subcluster
        # MGCs report edge status per subcluster via edge_subclusters.
        on_edge = all_subclusters[0].on_edge if len(all_subclusters) == 1 else False

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

        # Origin-spanning MGC on complete circular replicon(s): report the shortest
        # covering arc (start > end) when boundary subclusters actually wrap across
        # the origin within max_distance_bp (not two distant independent islands).
        if circular_contigs:
            if len(contigs_used) == 1 and circular_contigs.get(contigs_used[0]):
                final_mgc.is_circular = True
                ctg_len = (contig_lengths or {}).get(contigs_used[0], 0)
                wrapped = _origin_wrapped_span(
                    all_subclusters,
                    ctg_len,
                    max_origin_gap=self.max_distance_bp,
                    force_origin_bridge=(contigs_used[0] in origin_bridged_contigs),
                )
                if wrapped:
                    w_start, w_end, w_core_start, w_core_end, is_compact = wrapped
                    final_mgc.wraps_origin = True
                    final_mgc.wraps_origin_compact = is_compact
                    final_mgc.contig_length = ctg_len
                    final_mgc.origin_contigs[contigs_used[0]] = ctg_len
                    if is_compact:
                        final_mgc.start_bp = w_start
                        final_mgc.end_bp = w_end
                        final_mgc.cluster_core_start = w_core_start
                        final_mgc.cluster_core_end = w_core_end
            elif len(contigs_used) > 1:
                subs_by_ctg: Dict[str, List[MGCSubcluster]] = {}
                for sub in all_subclusters:
                    subs_by_ctg.setdefault(sub.contig_id, []).append(sub)
                for cid in contigs_used:
                    if not circular_contigs.get(cid):
                        continue
                    ctg_len = (contig_lengths or {}).get(cid, 0)
                    wrapped = _origin_wrapped_span(
                        subs_by_ctg.get(cid, []),
                        ctg_len,
                        max_origin_gap=self.max_distance_bp,
                        force_origin_bridge=(cid in origin_bridged_contigs),
                    )
                    if wrapped:
                        final_mgc.wraps_origin = True
                        final_mgc.origin_contigs[cid] = ctg_len
        return [final_mgc]
