"""High-performance profile-HMM homology search using PyHMMER."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pyhmmer
from ..data import get_hmm_path, load_rules
from .predict import PredictedGene

#: Bit-score cutoff used for a profile that has no ``trusted_cutoff`` in rules.yaml.
DEFAULT_SCORE_CUTOFF = 25.0
#: Default minimum fraction of the HMM profile covered by the alignment (0-1).
DEFAULT_HMM_COV_CUTOFF = 0.5


def _merged_span_length(intervals: List[Tuple[int, int]]) -> int:
    """Return total covered length of 1-based closed intervals after merging overlaps."""
    valid = sorted((min(a, b), max(a, b)) for a, b in intervals if a > 0 and b > 0)
    if not valid:
        return 0
    cur_s, cur_e = valid[0]
    total = 0
    for s, e in valid[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
    total += cur_e - cur_s + 1
    return total


def _compute_hit_coverages(hit, hmm_len: int, seq_len: int) -> Tuple[float, float]:
    """Compute ``(hmm_cov, seq_cov)`` in ``[0.0, 1.0]`` across aligned domain(s) of *hit*."""
    domains = [d for d in hit.domains if getattr(d, "included", True) and d.alignment is not None]
    if not domains:
        domains = [d for d in hit.domains if d.alignment is not None]
    if not domains and getattr(hit, "best_domain", None) is not None and hit.best_domain.alignment is not None:
        domains = [hit.best_domain]
    if not domains:
        return 0.0, 0.0

    hmm_intervals: List[Tuple[int, int]] = []
    seq_intervals: List[Tuple[int, int]] = []
    eff_hmm_len = hmm_len
    eff_seq_len = seq_len

    for d in domains:
        aln = d.alignment
        hmm_intervals.append((int(aln.hmm_from), int(aln.hmm_to)))
        seq_intervals.append((int(aln.target_from), int(aln.target_to)))
        if not eff_hmm_len and getattr(aln, "hmm_length", 0):
            eff_hmm_len = int(aln.hmm_length)
        if not eff_seq_len and getattr(aln, "target_length", 0):
            eff_seq_len = int(aln.target_length)

    hmm_cov = min(1.0, _merged_span_length(hmm_intervals) / eff_hmm_len) if eff_hmm_len > 0 else 0.0
    seq_cov = min(1.0, _merged_span_length(seq_intervals) / eff_seq_len) if eff_seq_len > 0 else 0.0
    return hmm_cov, seq_cov


@dataclass
class MagnetosomeHit:
    """Represents a significant HMM match for a gene."""
    gene_id: str
    contig_id: str
    profile: str
    parent_family: str
    score: float
    evalue: float
    role: str               # 'core', 'accessory', 'iron_context'
    lineage: str            # 'Alpha_classical', 'Desulfo_Mad', etc.
    color: str
    start: int
    end: int
    strand: int
    gene_index: int
    translation: str
    hmm_cov: float = 1.0
    seq_cov: float = 1.0

class MagnetosomeSearcher:
    """Orchestrates PyHMMER searches against the Magnetosome HMM library."""

    def __init__(
        self,
        hmm_path: Optional[Path] = None,
        rules_path: Optional[Path] = None,
        evalue_cutoff: float = 1e-5,
        score_cutoff: Optional[float] = None,
        hmm_cov_cutoff: float = DEFAULT_HMM_COV_CUTOFF,
        cpus: int = 0
    ):
        """
        Args:
            evalue_cutoff: maximum E-value of a reported hit.
            score_cutoff: global minimum bit score.  ``None`` (default) uses each
                profile's ``trusted_cutoff`` from rules.yaml (25.0 for all built-in
                profiles); an explicit value replaces the per-profile cutoffs for
                every profile.
            hmm_cov_cutoff: minimum fraction (0-1) of the HMM model covered by the hit.
            cpus: worker threads for PyHMMER (0 = all available cores).
        """
        if not (0.0 <= hmm_cov_cutoff <= 1.0):
            raise ValueError("hmm_cov_cutoff must be between 0.0 and 1.0")
        self.hmm_path = hmm_path or get_hmm_path()
        self.rules = load_rules(rules_path)
        self.evalue_cutoff = evalue_cutoff
        self.score_cutoff = score_cutoff
        self.hmm_cov_cutoff = float(hmm_cov_cutoff)
        self.cpus = cpus

        # Load HMMs into memory
        self.alphabet = pyhmmer.easel.Alphabet.amino()
        self.hmms = []
        with pyhmmer.plan7.HMMFile(str(self.hmm_path)) as hmm_file:
            self.hmms = list(hmm_file)

        # Profile metadata lookup from rules
        self.profile_meta = self.rules.get("profiles", {})

    def min_score_for(self, profile_name: str) -> float:
        """Effective bit-score cutoff for *profile_name*.

        An explicit ``score_cutoff`` (``-s/--score``) always wins; otherwise the
        profile's ``trusted_cutoff`` from rules.yaml is used.
        """
        if self.score_cutoff is not None:
            return float(self.score_cutoff)
        meta = self.profile_meta.get(profile_name, {})
        return float(meta.get("trusted_cutoff", DEFAULT_SCORE_CUTOFF))

    def search_genes(self, genes: List[PredictedGene]) -> Dict[str, MagnetosomeHit]:
        """Search a list of genes against magnetosome HMMs, returning best hit per gene."""
        if not genes or not self.hmms:
            return {}

        # Convert genes to digitized easel sequences
        easel_seqs = []
        gene_lookup = {}
        for g in genes:
            gene_lookup[g.gene_id] = g
            text_seq = pyhmmer.easel.TextSequence(
                name=g.gene_id.encode(),
                sequence=g.translation
            )
            easel_seqs.append(text_seq.digitize(self.alphabet))

        # Run hmmsearch
        # pyhmmer.hmmer.hmmsearch(queries, targets, cpus=cpus)
        search_results = list(
            pyhmmer.hmmer.hmmsearch(
                self.hmms,
                easel_seqs,
                cpus=self.cpus,
                E=self.evalue_cutoff
            )
        )

        # Process hits: for each gene, track (1) highest bit-score across all profiles
        # (to identify the best parent_family) and (2) best coverage-passing hit within
        # each parent_family (so domain sub-profiles like MamO_domain can win when
        # MamO_full has hmm_cov < hmm_cov_cutoff, without falling back to an unrelated
        # lower-scoring family).
        best_family_by_gene: Dict[str, Tuple[float, str]] = {}
        best_passing_by_gene_fam: Dict[Tuple[str, str], MagnetosomeHit] = {}

        for hmm, top_hits in zip(self.hmms, search_results):
            prof_name = hmm.name.decode() if isinstance(hmm.name, bytes) else str(hmm.name)
            meta = self.profile_meta.get(prof_name, {})
            parent_fam = meta.get("parent_family", prof_name)
            role = meta.get("role", "accessory")
            lineage = meta.get("lineage", "Other")
            color = meta.get("color", "#A8DADC")
            min_score = self.min_score_for(prof_name)

            raw_hmm_len = int(getattr(hmm, "M", 0))
            med_len = int(meta.get("median_length", 0) or 0)
            hmm_len = min(raw_hmm_len, med_len) if (raw_hmm_len > 0 and med_len > 0) else (raw_hmm_len or med_len)
            for hit in top_hits:
                g_id = hit.name.decode() if isinstance(hit.name, bytes) else str(hit.name)
                score = float(hit.score)
                evalue = float(hit.evalue)

                if evalue > self.evalue_cutoff or score < min_score:
                    continue

                if g_id not in best_family_by_gene or score > best_family_by_gene[g_id][0]:
                    best_family_by_gene[g_id] = (score, parent_fam)

                g_info = gene_lookup[g_id]
                hmm_cov, seq_cov = _compute_hit_coverages(hit, hmm_len, len(g_info.translation))
                if hmm_cov < self.hmm_cov_cutoff:
                    continue

                fam_key = (g_id, parent_fam)
                if fam_key in best_passing_by_gene_fam and score <= best_passing_by_gene_fam[fam_key].score:
                    continue

                best_passing_by_gene_fam[fam_key] = MagnetosomeHit(
                    gene_id=g_id,
                    contig_id=g_info.contig_id,
                    profile=prof_name,
                    parent_family=parent_fam,
                    score=score,
                    evalue=evalue,
                    role=role,
                    lineage=lineage,
                    color=color,
                    start=g_info.start,
                    end=g_info.end,
                    strand=g_info.strand,
                    gene_index=g_info.gene_index,
                    translation=g_info.translation,
                    hmm_cov=hmm_cov,
                    seq_cov=seq_cov,
                )

        best_hits_by_gene: Dict[str, MagnetosomeHit] = {}
        for g_id, (_, winning_fam) in best_family_by_gene.items():
            fam_key = (g_id, winning_fam)
            if fam_key in best_passing_by_gene_fam:
                best_hits_by_gene[g_id] = best_passing_by_gene_fam[fam_key]

        return best_hits_by_gene
