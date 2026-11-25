"""Fast, native prokaryotic gene prediction using Pyrodigal (Cython Prodigal)."""

from dataclasses import dataclass
from typing import List, Optional
from Bio.SeqRecord import SeqRecord
import pyrodigal

@dataclass
class PredictedGene:
    """Represents a gene predicted by Pyrodigal or parsed from existing annotation."""
    gene_id: str
    contig_id: str
    start: int          # 1-based start
    end: int            # 1-based end (inclusive)
    strand: int         # +1 or -1
    gene_index: int     # 0-based sequential order along the contig
    translation: str
    partial: bool = False
    score: float = 0.0
    rbs_motif: Optional[str] = None
    locus_tag: Optional[str] = None
    product: Optional[str] = None
    genome_name: Optional[str] = None

    @property
    def length_bp(self) -> int:
        return self.end - self.start + 1

    @property
    def length_aa(self) -> int:
        return len(self.translation)


class GenePredictor:
    """Pyrodigal gene prediction orchestrator."""

    def __init__(
        self,
        translation_table: int = 11,
        meta: Optional[bool] = None,
        min_gene_len: int = 90,
    ):
        self.translation_table = translation_table
        self.meta = meta
        self.min_gene_len = min_gene_len
        #: Mode actually used by the last ``predict_from_records`` call ('meta' or 'single').
        self.last_mode: Optional[str] = None
        #: Human-readable reason when single mode had to fall back to meta mode.
        self.last_fallback_reason: Optional[str] = None

    def predict_from_records(
        self,
        records: List[SeqRecord],
        genome_name: str = "genome"
    ) -> List[PredictedGene]:
        """Predict genes from a list of Biopython SeqRecords."""
        if not records:
            return []

        # Filter out records without defined nucleotide sequences (e.g. CONTIG-only GenBank)
        defined_records = []
        for r in records:
            try:
                if len(r.seq) > 0:
                    bytes(r.seq[:1])
                    defined_records.append(r)
            except Exception:
                continue
        if not defined_records:
            return []

        records = defined_records
        total_bp = sum(len(r.seq) for r in records)
        max_ctg_bp = max(len(r.seq) for r in records)

        # Determine single vs meta mode
        if self.meta is not None:
            use_meta = self.meta
        else:
            # Automatic: use meta mode for short/highly fragmented MAGs
            use_meta = total_bp < 100_000 or (len(records) > 200 and max_ctg_bp < 20_000)

        gf = pyrodigal.GeneFinder(meta=use_meta)
        self.last_fallback_reason = None

        # In single mode, train on the combined records
        if not use_meta:
            seq_bytes = [bytes(r.seq) for r in records if len(r.seq) >= 2000]
            if seq_bytes:
                try:
                    gf.train(*seq_bytes, translation_table=self.translation_table)
                except ValueError as e:
                    # Pyrodigal refuses to train on < 20 kb of sequence: fall back to meta mode
                    gf = pyrodigal.GeneFinder(meta=True)
                    use_meta = True
                    self.last_fallback_reason = (
                        f"single-genome training failed ({e}); used metagenome mode instead"
                    )
            else:
                # Fallback to meta if training sequences are too short
                gf = pyrodigal.GeneFinder(meta=True)
                use_meta = True
                self.last_fallback_reason = (
                    "no contig >= 2,000 bp available for single-genome training; used metagenome mode instead"
                )
        self.last_mode = "meta" if use_meta else "single"

        predicted_genes: List[PredictedGene] = []
        global_gene_idx = 1

        for rec in records:
            ctg_id = rec.id
            seq_bytes = bytes(rec.seq)
            if len(seq_bytes) < 30:
                continue

            genes = gf.find_genes(seq_bytes)
            ctg_gene_idx = 0
            for g in genes:
                trans = g.translate() if use_meta else g.translate(translation_table=self.translation_table)
                # Remove trailing stop if present
                trans = trans.rstrip("*")

                if len(trans) * 3 < self.min_gene_len:
                    continue

                gene_id = f"gene_{global_gene_idx:05d}"
                partial = (g.partial_begin != 0 or g.partial_end != 0)

                pg = PredictedGene(
                    gene_id=gene_id,
                    contig_id=ctg_id,
                    start=g.begin,
                    end=g.end,
                    strand=g.strand,
                    gene_index=ctg_gene_idx,
                    translation=trans,
                    partial=partial,
                    score=float(g.score),
                    rbs_motif=g.rbs_motif,
                    locus_tag=gene_id
                )
                predicted_genes.append(pg)
                ctg_gene_idx += 1
                global_gene_idx += 1

        return predicted_genes
