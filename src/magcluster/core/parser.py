"""Unified genome and annotation parser (FASTA, GenBank/GBFF, protein FASTA).

GFF3 files are *not* supported as input: they carry gene coordinates but no
protein or contig sequences.  Supply the matching GenBank/GBFF, nucleotide
FASTA or Prodigal-style protein FASTA instead.
"""

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
FAA_EXTENSIONS = {".faa", ".pep"}
GFF_EXTENSIONS = {".gff", ".gff3"}

#: Extensions picked up when a directory is given as input (optionally ``.gz``).
GENOME_EXTENSIONS = GENBANK_EXTENSIONS | FASTA_EXTENSIONS | FAA_EXTENSIONS

#: All extensions stripped from a file name to derive the genome ID.
KNOWN_EXTENSIONS = GENOME_EXTENSIONS | GFF_EXTENSIONS

PRODIGAL_HDR_RE = re.compile(
    r"^(?P<id>\S+)\s+#\s+(?P<start>\d+)\s+#\s+(?P<end>\d+)\s+#\s+(?P<strand>[-+]?\d+)\s+#\s+(?P<attr>.*)$"
)

# Header keywords that mark a FASTA record as a complete (circular) replicon.
CIRCULAR_HDR_RE = re.compile(r"complete\s+(genome|sequence)|\bcircular\b|topology=circular", re.IGNORECASE)

# Characters kept in genome IDs; everything else (whitespace, shell/filesystem
# specials) is replaced by '_'.  \w is Unicode-aware, so e.g. Chinese is kept.
_UNSAFE_NAME_RE = re.compile(r"[^\w.\-]+")


class InputFormatError(ValueError):
    """Raised when an input file cannot be parsed as a supported genome format."""

    def __init__(self, msg: str):
        # Third-party messages (e.g. Biopython) can span several lines: keep one line
        super().__init__(" ".join(str(msg).split()))


def _short_reason(exc: Exception, limit: int = 160) -> str:
    """Condense a third-party parser error to a short, single-line reason."""
    msg = " ".join(str(exc).split())
    if "comments at the beginning" in msg:  # Biopython >= 1.85 on leading non-'>' text
        return "text found before the first '>' header line"
    if len(msg) > limit:
        msg = msg[: limit - 3].rstrip() + "..."
    return msg


@dataclass
class ParsedGenome:
    """Represents a parsed genome with its contig sequences and genes."""
    name: str
    file_path: Path
    file_type: str                         # 'genbank', 'fasta', 'faa'
    contigs: Dict[str, SeqRecord] = field(default_factory=dict)
    genes: List[PredictedGene] = field(default_factory=list)
    needs_prediction: bool = False
    has_coordinates: bool = True
    circular_contigs: Dict[str, bool] = field(default_factory=dict)


def _strip_gz(name: str) -> str:
    return name[:-3] if name.lower().endswith(".gz") else name


def known_extension(path: Path) -> str:
    """Return the recognised (lower-case) genome extension of *path*, ignoring ``.gz``; '' if none."""
    ext = Path(_strip_gz(path.name)).suffix.lower()
    return ext if ext in KNOWN_EXTENSIONS else ""


def is_genome_file(path: Path) -> bool:
    """True if *path* has a supported genome extension (optionally gzipped)."""
    ext = known_extension(path)
    return bool(ext) and ext in GENOME_EXTENSIONS


def genome_name_from_path(path: Path) -> str:
    """Genome ID from a file name: the full name minus ``.gz`` and one known extension.

    ``strainX.v1.fna`` -> ``strainX.v1``; ``GCF_000513295.1.gbff.gz`` -> ``GCF_000513295.1``;
    files without a known extension keep their full name.
    """
    name = _strip_gz(path.name)
    ext = known_extension(path)
    if ext and name.lower().endswith(ext):
        name = name[: -len(ext)]
    return name or path.name


def sanitize_genome_name(name: str) -> str:
    """Make a genome ID safe for file names, GenBank LOCUS lines and Clinker.

    Whitespace and characters other than letters (including non-ASCII), digits,
    '.', '-' and '_' are replaced by '_'.  Leading dots are removed so that output
    files never become hidden files.
    """
    safe = _UNSAFE_NAME_RE.sub("_", name).lstrip(".")
    return safe or "genome"


def _open_text(path: Path):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    return opener(path, "rt", encoding="utf-8", errors="ignore")


def _check_not_binary(path: Path):
    """Reject binary (non-text) files with a clear message."""
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    try:
        with opener(path, "rb") as f:
            chunk = f.read(8192)
    except (OSError, EOFError) as e:
        raise InputFormatError(f"cannot read file ({_short_reason(e)})") from None
    if b"\x00" in chunk:
        raise InputFormatError("binary file, not a FASTA/GenBank text file")
    return chunk


def detect_file_type(path: Path) -> str:
    """Detect whether file is GenBank, nucleotide FASTA, protein FASTA (FAA) or GFF."""
    ext = known_extension(path)

    if ext in GENBANK_EXTENSIONS:
        return "genbank"
    if ext in FAA_EXTENSIONS:
        return "faa"
    if ext in FASTA_EXTENSIONS:
        return "fasta"
    if ext in GFF_EXTENSIONS:
        return "gff"

    # Inspect first non-empty line if extension is missing or unknown
    try:
        with _open_text(path) as f:
            for line in f:
                first_line = line.strip()
                if not first_line:
                    continue
                if first_line.startswith("LOCUS"):
                    return "genbank"
                if first_line.startswith(">"):
                    return "fasta"
                if first_line.startswith("##gff"):
                    return "gff"
                break
    except (OSError, EOFError):
        pass
    return "unknown"

def parse_genbank(gbk_path: Path, genome_name: Optional[str] = None) -> ParsedGenome:
    """Parse existing GenBank/GBFF file and extract already annotated CDS features.
    
    Robust to Biopython >= 1.80 (uses int(feature.location.start) rather than .position).
    """
    prefix = genome_name or sanitize_genome_name(genome_name_from_path(gbk_path))
    _check_not_binary(gbk_path)

    contigs = {}
    circular: Dict[str, bool] = {}
    genes: List[PredictedGene] = []
    gene_idx = 1

    try:
        with _open_text(gbk_path) as f:
            records = list(SeqIO.parse(f, "genbank"))
    except Exception as e:  # Biopython raises ValueError/AttributeError on malformed input
        raise InputFormatError(f"not a valid GenBank/GBFF file ({_short_reason(e)})") from None
    if not records:
        raise InputFormatError("no GenBank records found (expected a 'LOCUS' line)")

    for rec in records:
        contigs[rec.id] = rec
        circular[rec.id] = str(rec.annotations.get("topology", "")).lower() == "circular"
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
        needs_prediction=len(genes) == 0,
        circular_contigs=circular,
    )

def parse_fasta(fasta_path: Path, genome_name: Optional[str] = None) -> ParsedGenome:
    """Parse nucleotide FASTA, populating contigs, flagged for gene prediction."""
    prefix = genome_name or sanitize_genome_name(genome_name_from_path(fasta_path))
    _check_not_binary(fasta_path)

    contigs = {}
    circular: Dict[str, bool] = {}
    try:
        with _open_text(fasta_path) as f:
            for rec in SeqIO.parse(f, "fasta"):
                contigs[rec.id] = rec
                circular[rec.id] = bool(CIRCULAR_HDR_RE.search(rec.description))
    except Exception as e:  # Biopython raises ValueError on text that is not FASTA
        raise InputFormatError(f"not a valid FASTA file ({_short_reason(e)})") from None
    if not contigs:
        raise InputFormatError("no FASTA records found (expected a '>' header line)")
    if not any(len(rec.seq) for rec in contigs.values()):
        raise InputFormatError("FASTA records contain no sequence")

    return ParsedGenome(
        name=prefix,
        file_path=fasta_path,
        file_type="fasta",
        contigs=contigs,
        genes=[],
        needs_prediction=True,
        circular_contigs=circular,
    )

def parse_faa(faa_path: Path, genome_name: Optional[str] = None) -> ParsedGenome:
    """Parse protein FASTA, extracting coordinates if in Prodigal header format."""
    prefix = genome_name or sanitize_genome_name(genome_name_from_path(faa_path))
    _check_not_binary(faa_path)

    genes: List[PredictedGene] = []
    contig_counts = {}
    any_prodigal = False

    try:
        with _open_text(faa_path) as f:
            records = list(SeqIO.parse(f, "fasta"))
    except Exception as e:
        raise InputFormatError(f"not a valid protein FASTA file ({_short_reason(e)})") from None
    if not records:
        raise InputFormatError("no FASTA records found (expected a '>' header line)")

    for i, rec in enumerate(records):
        m = PRODIGAL_HDR_RE.match(rec.description)
        if m:
            any_prodigal = True
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
            start = 0
            end = 0
            strand = 0
            ctg_id = "unplaced"

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
        needs_prediction=False,
        has_coordinates=any_prodigal
    )

def parse_input(path: Path, genome_name: Optional[str] = None) -> ParsedGenome:
    """Automatically detect and parse input file into ParsedGenome.

    Raises :class:`InputFormatError` with a user-readable reason for unsupported
    or malformed files.
    """
    ftype = detect_file_type(path)
    if ftype == "genbank":
        return parse_genbank(path, genome_name)
    if ftype == "faa":
        return parse_faa(path, genome_name)
    if ftype == "fasta":
        return parse_fasta(path, genome_name)
    if ftype == "gff":
        raise InputFormatError(
            "GFF3 input is not supported (GFF3 has no protein/contig sequences); "
            "provide the matching GenBank/GBFF, nucleotide FASTA or Prodigal protein FASTA instead"
        )
    _check_not_binary(path)
    raise InputFormatError(
        "unrecognised file format (expected GenBank/GBFF starting with 'LOCUS' or FASTA starting with '>')"
    )
