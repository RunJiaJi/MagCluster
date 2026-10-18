"""High-performance profile-HMM homology search using PyHMMER."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional
import pyhmmer
from ..data import get_hmm_path, load_rules
from .predict import PredictedGene

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

class MagnetosomeSearcher:
    """Orchestrates PyHMMER searches against the Magnetosome HMM library."""

    def __init__(
        self,
        hmm_path: Optional[Path] = None,
        rules_path: Optional[Path] = None,
        evalue_cutoff: float = 1e-5,
        score_cutoff: float = 25.0,
        cpus: int = 0
    ):
        self.hmm_path = hmm_path or get_hmm_path()
        self.rules = load_rules(rules_path)
        self.evalue_cutoff = evalue_cutoff
        self.score_cutoff = score_cutoff
        self.cpus = cpus

        # Load HMMs into memory
        self.alphabet = pyhmmer.easel.Alphabet.amino()
        self.hmms = []
        with pyhmmer.plan7.HMMFile(str(self.hmm_path)) as hmm_file:
            self.hmms = list(hmm_file)

        # Profile metadata lookup from rules
        self.profile_meta = self.rules.get("profiles", {})

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

        # Process hits: find best profile for each gene
        best_hits_by_gene: Dict[str, MagnetosomeHit] = {}

        for hmm, top_hits in zip(self.hmms, search_results):
            prof_name = hmm.name.decode() if isinstance(hmm.name, bytes) else str(hmm.name)
            meta = self.profile_meta.get(prof_name, {})
            parent_fam = meta.get("parent_family", prof_name)
            role = meta.get("role", "accessory")
            lineage = meta.get("lineage", "Other")
            color = meta.get("color", "#A8DADC")
            min_score = meta.get("trusted_cutoff", self.score_cutoff)

            for hit in top_hits:
                g_id = hit.name.decode() if isinstance(hit.name, bytes) else str(hit.name)
                score = float(hit.score)
                evalue = float(hit.evalue)

                if evalue > self.evalue_cutoff or score < min_score:
                    continue

                # If gene already hit a profile, keep the higher scoring one
                if g_id in best_hits_by_gene:
                    if score <= best_hits_by_gene[g_id].score:
                        continue

                g_info = gene_lookup[g_id]
                best_hits_by_gene[g_id] = MagnetosomeHit(
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
                    translation=g_info.translation
                )

        return best_hits_by_gene
