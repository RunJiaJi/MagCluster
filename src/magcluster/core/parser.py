"""Unified genome and annotation parser (FASTA, GenBank/GBFF, GFF3, FAA)."""

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple, Union
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from .predict import PredictedGene

GENBANK_EXTENSIONS = {".gbk", ".gbff", ".gbf", ".gb"}
FASTA_EXTENSIONS = {".fa", ".fasta", ".fna", ".fsa", ".fsa_nt"}
FAA_EXTENSIONS = {".faa", ".pep", ".faa.gz"}

PRODIGAL_HDR_RE = re.compile(
    r"^(?P<id>\S+)\s+#\s+(?P<start>\d+)\s+#\s+(?P<end>\d+)\s+#\s+(?P<strand>[-+]?\d+)\s+#\s+(?P<attr>.*)$"
)

@dataclass
class ParsedGenome:
    """Represents a parsed genome with its contig sequences and genes."""
    name: str
    file_path: Path
    file_type: str                         # 'genbank', 'fasta', 'faa', 'gff'
    contigs: Dict[str, SeqRecord] = field(default_factory=dict)
    genes: List[PredictedGene] = field(default_factory=list)
    needs_prediction: bool = False

def detect_file_type(path: Path) -> str:
    """Detect whether file is GenBank, FASTA nucleotide, or FAA."""
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    ext = Path(name).suffix

    if ext in GENBANK_EXTENSIONS:
        return "genbank"
    if ext in FAA_EXTENSIONS:
        return "faa"
    if ext in FASTA_EXTENSIONS:
        return "fasta"
    if ext in {".gff", ".gff3"}:
        return "gff"

    # Inspect first line if extension is ambiguous
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="ignore") as f:
            first_line = f.readline().strip()
            if first_line.startswith("LOCUS"):
                return "genbank"
            if first_line.startswith(">"):
                return "fasta"
            if first_line.startswith("##gff"):
                return "gff"
    except Exception:
        pass
    return "fasta"

def parse_genbank(gbk_path: Path, genome_name: Optional[str] = None) -> ParsedGenome:
    """Parse existing GenBank/GBFF file and extract already annotated CDS features.
    
    Robust to Biopython >= 1.80 (uses int(feature.location.start) rather than .position).
    """
    prefix = genome_name or gbk_path.stem.split(".")[0]
    opener = gzip.open if str(gbk_path).endswith(".gz") else open

    contigs = {}
    genes: List[PredictedGene] = []
    gene_idx = 1

    with opener(gbk_path, "rt", encoding="utf-8", errors="ignore") as f:
        records = list(SeqIO.parse(f, "genbank"))

    for rec in records:
        contigs[rec.id] = rec
        ctg_gene_idx = 0

        for feat in rec.features:
            if feat.type != "CDS":
                continue

            # Safely extract coordinates (compatible with all Biopython versions)
            start_pos = int(feat.location.start) + 1  # convert to 1-based closed
            end_pos = int(feat.location.end)          # 1-based closed
            strand = 1 if feat.location.strand >= 0 else -1

            # Extract translation
            trans = feat.qualifiers.get("translation", [None])[0]
            if not trans:
                # If translation is missing, attempt in silico translation if nucleotide sequence exists
                if len(rec.seq) > 0:
                    try:
                        sub_seq = feat.extract(rec.seq)
                        trans = str(sub_seq.translate()).rstrip("*")
                    except Exception:
                        continue
                else:
                    continue

            locus_tag = feat.qualifiers.get("locus_tag", [f"{rec.id}_{gene_idx}"])[0]
            product = feat.qualifiers.get("product", [None])[0]
            gene_id = locus_tag

            pg = PredictedGene(
                gene_id=gene_id,
                contig_id=rec.id,
                start=start_pos,
                end=end_pos,
                strand=strand,
                gene_index=ctg_gene_idx,
                translation=trans,
                locus_tag=locus_tag,
                product=product
            )
            genes.append(pg)
            ctg_gene_idx += 1
            gene_idx += 1

    return ParsedGenome(
        name=prefix,
        file_path=gbk_path,
        file_type="genbank",
        contigs=contigs,
        genes=genes,
        needs_prediction=len(genes) == 0
    )

def parse_fasta(fasta_path: Path, genome_name: Optional[str] = None) -> ParsedGenome:
    """Parse nucleotide FASTA, populating contigs, flagged for gene prediction."""
    prefix = genome_name or fasta_path.stem.split(".")[0]
    opener = gzip.open if str(fasta_path).endswith(".gz") else open

    contigs = {}
    with opener(fasta_path, "rt", encoding="utf-8", errors="ignore") as f:
        for rec in SeqIO.parse(f, "fasta"):
            contigs[rec.id] = rec

    return ParsedGenome(
        name=prefix,
        file_path=fasta_path,
        file_type="fasta",
        contigs=contigs,
        genes=[],
        needs_prediction=True
    )

def parse_faa(faa_path: Path, genome_name: Optional[str] = None) -> ParsedGenome:
    """Parse protein FASTA, extracting coordinates if in Prodigal header format."""
    prefix = genome_name or faa_path.stem.split(".")[0]
    opener = gzip.open if str(faa_path).endswith(".gz") else open

    genes: List[PredictedGene] = []
    contig_counts = {}

    with opener(faa_path, "rt", encoding="utf-8", errors="ignore") as f:
        for i, rec in enumerate(SeqIO.parse(f, "fasta")):
            m = PRODIGAL_HDR_RE.match(rec.description)
            if m:
                gene_id = m["id"]
                start = int(m["start"])
                end = int(m["end"])
                strand = int(m["strand"])
                attr = m["attr"]
                # Try to extract contig from attributes or header
                ctg_match = re.search(r"contig=([^;]+)", attr)
                ctg_id = ctg_match.group(1) if ctg_match else gene_id.rsplit("_", 1)[0]
            else:
                gene_id = rec.id
                start = 1 + i * 1000
                end = start + len(rec.seq) * 3
                strand = 1
                ctg_id = "contig_1"

            c_idx = contig_counts.get(ctg_id, 0)
            contig_counts[ctg_id] = c_idx + 1

            pg = PredictedGene(
                gene_id=gene_id,
                contig_id=ctg_id,
                start=start,
                end=end,
                strand=strand,
                gene_index=c_idx,
                translation=str(rec.seq).rstrip("*"),
                locus_tag=gene_id
            )
            genes.append(pg)

    return ParsedGenome(
        name=prefix,
        file_path=faa_path,
        file_type="faa",
        contigs={},
        genes=genes,
        needs_prediction=False
    )

def parse_input(path: Path, genome_name: Optional[str] = None) -> ParsedGenome:
    """Automatically detect and parse input file into ParsedGenome."""
    ftype = detect_file_type(path)
    if ftype == "genbank":
        return parse_genbank(path, genome_name)
    if ftype == "faa":
        return parse_faa(path, genome_name)
    # Default to nucleotide fasta
    return parse_fasta(path, genome_name)
