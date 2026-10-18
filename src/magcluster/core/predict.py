"""Fast, native prokaryotic gene prediction using Pyrodigal (Cython Prodigal)."""

import gzip
import io
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple, Union
from Bio import SeqIO
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

    @property
    def length_bp(self) -> int:
        return self.end - self.start + 1

    @property
    def length_aa(self) -> int:
        return len(self.translation)

    def to_gff3_line(self, source: str = "pyrodigal") -> str:
        strand_char = "+" if self.strand >= 0 else "-"
        attrs = f"ID={self.gene_id};partial={'true' if self.partial else 'false'}"
        if self.locus_tag:
            attrs += f";locus_tag={self.locus_tag}"
        if self.product:
            attrs += f";product={self.product}"
        return f"{self.contig_id}\t{source}\tCDS\t{self.start}\t{self.end}\t{self.score:.1f}\t{strand_char}\t0\t{attrs}"

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

    def _open_fasta(self, fasta_path: Path) -> Iterator[SeqRecord]:
        """Open plain or gzipped FASTA."""
        if str(fasta_path).endswith(".gz"):
            with gzip.open(fasta_path, "rt") as f:
                yield from SeqIO.parse(f, "fasta")
        else:
            with open(fasta_path, "rt") as f:
                yield from SeqIO.parse(f, "fasta")

    def predict_from_records(
        self,
        records: List[SeqRecord],
        genome_name: str = "genome"
    ) -> List[PredictedGene]:
        """Predict genes from a list of Biopython SeqRecords."""
        if not records:
            return []

        total_bp = sum(len(r.seq) for r in records)
        max_ctg_bp = max(len(r.seq) for r in records)

        # Determine single vs meta mode
        if self.meta is not None:
            use_meta = self.meta
        else:
            # Automatic: use meta mode for short/highly fragmented MAGs
            use_meta = total_bp < 100_000 or (len(records) > 200 and max_ctg_bp < 20_000)

        gf = pyrodigal.GeneFinder(meta=use_meta)
        
        # In single mode, train on the combined records
        if not use_meta:
            seq_bytes = [bytes(r.seq) for r in records if len(r.seq) >= 2000]
            if seq_bytes:
                gf.train(*seq_bytes, translation_table=self.translation_table)
            else:
                # Fallback to meta if training sequences are too short
                gf = pyrodigal.GeneFinder(meta=True)

        predicted_genes: List[PredictedGene] = []
        global_gene_idx = 1

        for rec in records:
            ctg_id = rec.id
            seq_bytes = bytes(rec.seq)
            if len(seq_bytes) < 30:
                continue

            genes = gf.find_genes(seq_bytes)
            for ctg_gene_idx, g in enumerate(genes):
                trans = g.translate(translation_table=self.translation_table)
                # Remove trailing stop if present
                trans = trans.rstrip("*")

                if len(trans) * 3 < self.min_gene_len:
                    continue

                gene_id = f"{genome_name}_{global_gene_idx:05d}"
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
                global_gene_idx += 1

        return predicted_genes

    def predict_from_file(
        self,
        fasta_path: Path,
        genome_name: Optional[str] = None
    ) -> List[PredictedGene]:
        """Predict genes directly from a FASTA file."""
        records = list(self._open_fasta(fasta_path))
        prefix = genome_name or fasta_path.name.split(".")[0]
        return self.predict_from_records(records, genome_name=prefix)

    @staticmethod
    def write_faa(genes: List[PredictedGene], out_path: Path):
        """Export predicted proteins to FASTA format."""
        with open(out_path, "w", encoding="utf-8") as f:
            for g in genes:
                hdr = f">{g.gene_id} # {g.start} # {g.end} # {g.strand} # ID={g.gene_id};contig={g.contig_id};idx={g.gene_index}"
                f.write(f"{hdr}\n{g.translation}\n")

    @staticmethod
    def write_gff3(genes: List[PredictedGene], out_path: Path):
        """Export predicted genes to GFF3 format."""
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("##gff-version 3\n")
            for g in genes:
                f.write(f"{g.to_gff3_line()}\n")
