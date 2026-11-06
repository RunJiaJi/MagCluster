"""Clinker-compatible GenBank slice exporter and automated visualization runner."""

import gzip
import json
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from ..core.cluster import MGCCluster
from ..core.parser import InputFormatError, _check_not_binary, _short_reason, parse_genbank
from ..core.predict import PredictedGene
from ..core.search import MagnetosomeHit, MagnetosomeSearcher
from ..data import load_rules
from ..utils.logger import log

#: Maximum length of ACCESSION / VERSION content in Biopython InsdcIO (80 - 12 = 68).
MAX_GBK_ID_LEN = 68
#: Standard GenBank LOCUS identifier column width.
MAX_LOCUS_LEN = 16
#: Default maximum number of MGC clusters plotted per Clinker HTML file.
DEFAULT_CLINKER_BATCH_SIZE = 50

# Standard INSDC / RefSeq assembly or sequence accession embedded in genome names.
_ACCESSION_RE = re.compile(
    r"(?:^|(?<=_))"
    r"(GC[AF]_\d+(?:\.\d+)?"
    r"|(?:NC|NZ|NW|NT|NG|CP|AP|AE|AL|AM|BA|BK|BX|CR|CT|CU|FM|FN|FO|FP|FQ|FR|HE|HF|HG|LN|LO|LP|LQ|LR|LS|LT|OU|OV|OW|OX|OY|OZ)_?[A-Z0-9]+(?:\.\d+)?"
    r"|[A-Z]{4,6}\d{8,}(?:\.\d+)?)"
    r"(?=$|[_.-])"
)

# Metadata written into CDS /note qualifiers by export_cluster_to_genbank().
_MAGCLUSTER_NOTE_RE = re.compile(
    r"Profile=([^;\s]+)\s*;\s*Family=([^;\s]+)"
    r"(?:\s*;\s*Score=([0-9.eE+-]+))?"
    r"(?:\s*;\s*Evalue=([0-9.eE+-]+))?"
)

_CACHED_PROFILE_META: Optional[Dict[str, dict]] = None
_CACHED_SEARCHERS: Dict[int, MagnetosomeSearcher] = {}


def _profile_meta_map() -> Dict[str, dict]:
    global _CACHED_PROFILE_META
    if _CACHED_PROFILE_META is None:
        _CACHED_PROFILE_META = load_rules().get("profiles", {})
    return _CACHED_PROFILE_META


def _get_cached_searcher(cpus: int) -> MagnetosomeSearcher:
    if cpus not in _CACHED_SEARCHERS:
        _CACHED_SEARCHERS[cpus] = MagnetosomeSearcher(cpus=cpus)
    return _CACHED_SEARCHERS[cpus]

# Classic MagCluster v0.2.12 color palette extended for novel families
FAMILY_COLOR_PALETTE = {
    "mam": "#EA4335",   # Classic MagCluster Red (v0.2.12)
    "mms": "#FBBC05",   # Classic MagCluster Yellow (v0.2.12)
    "mad": "#34A853",   # Classic MagCluster Green (v0.2.12)
    "man": "#4385F4",   # Classic MagCluster Blue (v0.2.12)
    "feo": "#FFE066",   # Iron metabolism / transport (v0.2.12)
    "mae": "#8338EC",   # Vibrant Purple
    "maq": "#06D6A0",   # Vibrant Teal
    "mca": "#FF006E",   # Vibrant Magenta / Pink
    "Other": "#BDBDBD",  # Muted grey for flanking/intervening non-MGC genes
}

def get_family_class(gene_name: str, hit: Optional[MagnetosomeHit] = None) -> str:
    """Classify magnetosome gene into parent prefix group (mam, mms, mad, man, mae, maq, mca, feo, Other)."""
    if hit:
        if hit.role == "iron_context" or hit.parent_family.startswith("Feo"):
            return "feo"
        fam_lower = hit.parent_family.lower()
        for prefix in ["mam", "mms", "mad", "man", "mae", "maq", "mca"]:
            if fam_lower.startswith(prefix):
                return prefix
        return "Other"

    gn = gene_name.lower()
    for prefix in ["mam", "mms", "mad", "man", "mae", "maq", "mca", "feo"]:
        if prefix in gn:
            return prefix
    if "iron" in gn:
        return "feo"
    return "Other"

def get_display_label(hit: MagnetosomeHit) -> str:
    """Return the display label for a magnetosome hit in figures.
    
    Retains specific variant/domain labels (e.g. MamE_Cter, MamQ_greigite, MamO_domain),
    while omitting the '_full' suffix for full-length proteins (e.g. MamE_full -> MamE).
    """
    if hit.profile.endswith("_full"):
        return hit.profile[:-5]
    return hit.profile


def shorten_gbk_record_ids(
    genome_name: str,
    cluster_id: str,
    slice_idx: int = 1,
    n_slices: int = 1,
) -> Tuple[str, str]:
    """Return ``(rec_id, rec_name)`` for a sliced GenBank record.

    When the full identifier fits within Biopython's 68-character ``ACCESSION`` /
    ``VERSION`` line limit, the original v1.0 identifier and ``[:16]`` LOCUS name
    are kept unchanged.  When the identifier exceeds 68 characters, an assembly or
    sequence accession (e.g. ``GCA_002753255.1``) replaces the full genome name—
    falling back to safe truncation if no standard accession is present—so Biopython
    never emits header length warnings while the ``.gbk`` file name stays unchanged.
    """
    raw_id = f"{genome_name}_loc{slice_idx}" if n_slices > 1 else cluster_id
    raw_id = re.sub(r"\s+", "_", raw_id)
    if len(raw_id) <= MAX_GBK_ID_LEN:
        return raw_id, raw_id[:MAX_LOCUS_LEN]

    suffix = f"_loc{slice_idx}" if n_slices > 1 else "_MGC"
    clean_genome = re.sub(r"\s+", "_", genome_name)
    m = _ACCESSION_RE.search(clean_genome)
    base = m.group(1) if m else clean_genome

    max_base = max(1, MAX_GBK_ID_LEN - len(suffix))
    if len(base) > max_base:
        base = base[:max_base].rstrip("_.-") or base[:max_base]

    rec_id = f"{base}{suffix}"
    if len(rec_id) <= MAX_LOCUS_LEN:
        rec_name = rec_id
    elif n_slices > 1 and len(suffix) < MAX_LOCUS_LEN:
        rec_name = f"{base[:MAX_LOCUS_LEN - len(suffix)]}{suffix}"
    else:
        rec_name = base[:MAX_LOCUS_LEN]
    return rec_id, rec_name


def _wrap_genbank_definition(text: str, max_token: int = MAX_GBK_ID_LEN) -> str:
    """Ensure no single whitespace-delimited token in DEFINITION exceeds *max_token* chars.

    Biopython's ``GenBankWriter._write_multi_line`` wraps ``DEFINITION`` only on spaces
    and emits ``BiopythonWarning: Annotation '...' too long`` if any single token exceeds
    68 characters (80 - 12). Short descriptions are returned byte-for-byte unchanged.
    """
    out_words: List[str] = []
    for word in text.split(" "):
        if len(word) <= max_token:
            out_words.append(word)
            continue
        for part in word.replace("_", " ").split(" "):
            while len(part) > max_token:
                out_words.append(part[:max_token])
                part = part[max_token:]
            if part:
                out_words.append(part)
    return " ".join(out_words)


def _hit_from_profile_note(
    note_text: str,
    gene_id: str = "",
    contig_id: str = "",
    start: int = 0,
    end: int = 0,
    strand: int = 1,
    gene_index: int = 0,
    translation: str = "",
) -> Optional[MagnetosomeHit]:
    """Reconstruct a :class:`MagnetosomeHit` from a MagCluster CDS ``/note`` string."""
    if not note_text:
        return None
    m = _MAGCLUSTER_NOTE_RE.search(note_text)
    if not m:
        return None
    profile = m.group(1)
    parent_fam = m.group(2)
    score = float(m.group(3)) if m.group(3) else 0.0
    evalue = float(m.group(4)) if m.group(4) else 0.0
    meta = _profile_meta_map().get(profile, {})
    default_role = "iron_context" if parent_fam.startswith("Feo") else "core"
    return MagnetosomeHit(
        gene_id=gene_id,
        contig_id=contig_id,
        profile=profile,
        parent_family=parent_fam,
        score=score,
        evalue=evalue,
        role=meta.get("role", default_role),
        lineage=meta.get("lineage", "Other"),
        color=meta.get("color", "#EA4335"),
        start=start,
        end=end,
        strand=strand,
        gene_index=gene_index,
        translation=translation,
    )


def extract_hits_from_genbank(
    gbk_path: Path,
    searcher: Optional[MagnetosomeSearcher] = None,
    cpus: int = 0,
) -> Dict[str, MagnetosomeHit]:
    """Extract magnetosome hit annotations from a sliced or external GenBank file.

    For GenBank files produced by MagCluster (``clusters_gbk/*.gbk``), hit metadata
    (``Profile``, ``Family``, ``Score``, ``Evalue``) is read directly from CDS
    ``/note`` qualifiers without re-running HMM searches.  For external GenBank
    files that carry no MagCluster ``/note`` metadata, CDS translations are searched
    against the built-in 92 magnetosome HMMs using *searcher* (lazily initialized if
    not supplied).
    """
    _check_not_binary(gbk_path)
    opener = gzip.open if str(gbk_path).lower().endswith(".gz") else open
    try:
        with opener(gbk_path, "rt", encoding="utf-8", errors="ignore") as f:
            records = list(SeqIO.parse(f, "genbank"))
    except Exception as e:
        raise InputFormatError(f"not a valid GenBank/GBFF file ({_short_reason(e)})") from None
    if not records:
        raise InputFormatError("no GenBank records found (expected a 'LOCUS' line)")

    is_magcluster_slice = any("Magnetosome Gene Cluster" in (r.description or "") for r in records)
    note_hits: Dict[str, MagnetosomeHit] = {}

    gene_idx = 0
    for rec in records:
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            gene_idx += 1
            lt = feat.qualifiers.get("locus_tag", [None])[0]
            gn = feat.qualifiers.get("gene", [None])[0]
            key = lt or gn or f"{rec.id}_{gene_idx}"
            notes = " ".join(feat.qualifiers.get("note", []))
            h = _hit_from_profile_note(
                notes,
                gene_id=key,
                contig_id=rec.id,
                start=int(feat.location.start) + 1,
                end=int(feat.location.end),
                strand=1 if (feat.location.strand or 1) >= 0 else -1,
                gene_index=gene_idx - 1,
                translation=feat.qualifiers.get("translation", [""])[0],
            )
            if h is not None:
                note_hits[key] = h

    if note_hits or is_magcluster_slice:
        return note_hits

    # External GenBank file: annotate CDS features with the built-in magnetosome HMMs
    parsed = parse_genbank(gbk_path)
    if not parsed.genes:
        return {}
    active_searcher = searcher if searcher is not None else _get_cached_searcher(cpus)
    return active_searcher.search_genes(parsed.genes)


def export_cluster_to_genbank(
    cluster: MGCCluster,
    contigs_dict: Dict[str, SeqRecord],
    hits: Dict[str, MagnetosomeHit],
    out_path: Path
):
    """Export sliced GenBank file for an MGC cluster (supporting discrete subclusters & contigs)."""
    subclusters = getattr(cluster, "subclusters", None) or []
    if subclusters:
        slices = [
            (sub.contig_id, sorted(sub.all_genes, key=lambda x: x.start))
            for sub in subclusters
        ]
    else:
        # Fallback grouping by contig if subclusters not populated
        genes_by_ctg: Dict[str, List[PredictedGene]] = {}
        for g in cluster.all_genes:
            genes_by_ctg.setdefault(g.contig_id, []).append(g)
        slices = [
            (ctg_id, sorted(ctg_genes, key=lambda x: x.start))
            for ctg_id, ctg_genes in genes_by_ctg.items()
        ]

    records_to_write = []

    for i, (ctg_id, ctg_genes) in enumerate(slices, 1):
        if not ctg_genes:
            continue

        offset = ctg_genes[0].start - 1  # 0-based slice start
        end_idx = ctg_genes[-1].end      # 0-based slice end (exclusive)

        ctg_rec = contigs_dict.get(ctg_id)
        if ctg_rec and len(ctg_rec.seq) >= end_idx:
            sub_seq = ctg_rec.seq[offset:end_idx]
        else:
            span_len = max(1, end_idx - offset)
            sub_seq = Seq("N" * span_len)

        rec_id, rec_name = shorten_gbk_record_ids(
            cluster.genome_name, cluster.cluster_id, slice_idx=i, n_slices=len(slices)
        )
        desc_label = f"slice {i}" if len(slices) > 1 else "cluster"
        raw_desc = (
            f"Magnetosome Gene Cluster {desc_label} from {cluster.genome_name} "
            f"({ctg_id}:{offset+1}-{end_idx})"
        )
        new_rec = SeqRecord(
            sub_seq,
            id=rec_id,
            name=rec_name,
            description=_wrap_genbank_definition(raw_desc),
        )
        new_rec.annotations["molecule_type"] = "DNA"
        new_rec.annotations["topology"] = "linear"

        features = []
        # Source feature
        features.append(SeqFeature(
            FeatureLocation(0, len(sub_seq)),
            type="source",
            qualifiers={"organism": [cluster.genome_name]}
        ))

        # CDS features
        for g in ctg_genes:
            feat_start = max(0, g.start - 1 - offset)
            feat_end = min(len(sub_seq), g.end - offset)
            feat_strand = g.strand
            if feat_end <= feat_start:
                continue

            loc = FeatureLocation(feat_start, feat_end, strand=feat_strand)
            qualifiers = {
                "locus_tag": [g.gene_id],
                "translation": [g.translation]
            }

            if g.gene_id in hits:
                h = hits[g.gene_id]
                protein_name = get_display_label(h)
                qualifiers["product"] = [protein_name]
                qualifiers["gene"] = [protein_name]
                qualifiers["label"] = [protein_name]
                qualifiers["protein"] = [protein_name]
                qualifiers["note"] = [f"Profile={h.profile}; Family={h.parent_family}; Score={h.score:.1f}; Evalue={h.evalue:.2e}"]
            else:
                qualifiers["product"] = [g.product or "hypothetical protein"]
                qualifiers["gene"] = [""]
                qualifiers["label"] = [""]

            features.append(SeqFeature(loc, type="CDS", qualifiers=qualifiers))

        new_rec.features = features
        records_to_write.append(new_rec)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records_to_write, str(out_path), "genbank")

def export_all_clusters_for_clinker(
    clusters: List[MGCCluster],
    contigs_lookup: Dict[str, Dict[str, SeqRecord]],
    hits: Dict[str, MagnetosomeHit],
    out_dir: Path,
    hits_by_genome: Optional[Dict[str, Dict[str, MagnetosomeHit]]] = None,
) -> List[Path]:
    """Export sliced GenBank files for clinker visualization without creating extra CSV tables.

    ``hits_by_genome`` (genome name -> hits of that genome) takes precedence over the
    flat ``hits`` dict so that identical gene IDs in different genomes never mix.
    """
    gbk_dir = out_dir / "clusters_gbk"
    gbk_dir.mkdir(parents=True, exist_ok=True)
    for old_gbk in gbk_dir.glob("*.gbk"):
        old_gbk.unlink()

    gbk_files = []
    for c in clusters:
        gbk_file = gbk_dir / f"{c.cluster_id}.gbk"
        genome_contigs = contigs_lookup.get(c.genome_name, {})
        g_hits = hits_by_genome.get(c.genome_name, {}) if hits_by_genome is not None else hits
        export_cluster_to_genbank(c, genome_contigs, g_hits, gbk_file)
        gbk_files.append(gbk_file)

    return gbk_files

def customize_clinker_html(
    plot_html: Path,
    hits: Dict[str, MagnetosomeHit],
    hits_by_cluster: Optional[Dict[str, Dict[str, MagnetosomeHit]]] = None,
) -> bool:
    """Post-process Clinker HTML to apply v0.2.12 family colors and protein name labels.

    ``hits_by_cluster`` maps a Clinker cluster name (the GenBank file stem, i.e. the
    MGC ID) to the hits of that genome; clusters not in the map use ``hits``.
    """
    if not plot_html.exists():
        return False

    with open(plot_html, "r", encoding="utf-8") as f:
        html = f.read()

    m = re.search(r"<script>const data=(.*?);function serialise", html)
    if not m:
        log.warning("Could not locate Clinker data JSON in HTML; skipped custom coloring.")
        return False

    try:
        data = json.loads(m.group(1))
    except Exception as e:
        log.warning(f"Failed to parse Clinker data JSON: {e}")
        return False

    # Initialize groups
    group_dict = {
        cls_name: {
            "uid": f"group_{cls_name}",
            "label": cls_name,
            "genes": [],
            "hidden": False,
            "colour": color
        }
        for cls_name, color in FAMILY_COLOR_PALETTE.items()
    }

    # Iterate over all genes in every cluster locus
    for cluster in data.get("clusters", []):
        c_hits = hits
        if hits_by_cluster is not None:
            c_hits = hits_by_cluster.get(cluster.get("name"), hits)
        by_coord = {
            (h.start, h.end, h.strand): h
            for h in c_hits.values()
            if h.end > 0
        }
        for locus in cluster.get("loci", []):
            for gene in locus.get("genes", []):
                names = gene.get("names", {})
                lt = names.get("locus_tag")
                gn = names.get("gene")

                hit = None
                if lt and lt in c_hits:
                    hit = c_hits[lt]
                elif gn and gn in c_hits:
                    hit = c_hits[gn]
                elif names.get("note"):
                    hit = _hit_from_profile_note(names.get("note", ""), gene_id=lt or gn or "")
                elif by_coord and "start" in gene and "end" in gene:
                    hit = by_coord.get((int(gene["start"]) + 1, int(gene["end"]), int(gene.get("strand", 1))))

                if hit:
                    protein_name = get_display_label(hit)
                    gene["label"] = protein_name
                    gene["is_mag"] = True
                    gene.setdefault("names", {})["protein"] = protein_name
                    gene["names"]["gene"] = protein_name
                    gene["names"]["product"] = protein_name
                    gene["names"]["label"] = protein_name

                    f_class = get_family_class(hit.parent_family, hit)
                    group_dict[f_class]["genes"].append(gene["uid"])
                else:
                    # Non-magnetosome gene: blank label to prevent arrow clutter
                    gene["label"] = ""
                    gene["is_mag"] = False
                    gene["names"] = {
                        "protein": "",
                        "gene": "",
                        "label": "",
                        "locus_tag": "",
                        "product": ""
                    }
                    group_dict["Other"]["genes"].append(gene["uid"])

    # Retain groups that have at least one gene
    data["groups"] = [
        g for g in group_dict.values()
        if len(g["genes"]) > 0
    ]

    new_html = html[:m.start(1)] + json.dumps(data) + html[m.end(1):]

    # Prevent clustermap.js from falling back to t.uid when t.label is empty
    new_html = new_html.replace("labelText:t=>t.label||t.uid", "labelText:t=>t.label||\"\"")

    # Prevent sidebar label selector from ever applying labels to non-magnetosome genes
    new_html = new_html.replace(
        ".each(d => d.label = d.names[type] ? d.names[type] : d.label)",
        ".each(d => d.label = d.is_mag ? (d.names[type] || d.label) : \"\")"
    )

    with open(plot_html, "w", encoding="utf-8") as f:
        f.write(new_html)

    return True

class ClinkerError(RuntimeError):
    """Raised when Clinker is missing or fails to produce the HTML figure."""


def parse_clinker_opts(clinker_opts: Optional[str]) -> List[str]:
    """Split the ``--clinker-opts`` string like a shell would.

    Raises ValueError with a readable message (e.g. unbalanced quotes).
    """
    if not clinker_opts:
        return []
    try:
        return shlex.split(clinker_opts)
    except ValueError as e:
        raise ValueError(f"cannot parse --clinker-opts {clinker_opts!r}: {e}") from None


def find_clinker() -> Optional[str]:
    """Locate the clinker executable, preferring the one next to the running interpreter."""
    env_bin = Path(sys.executable).parent / "clinker"
    if env_bin.exists():
        return str(env_bin)
    return shutil.which("clinker")


def run_clinker_visualization(
    gbk_files: List[Path],
    plot_html: Path,
    hits: Dict[str, MagnetosomeHit],
    jobs: int = 0,
    identity: float = 0.3,
    clinker_opts: Optional[str] = None,
    hits_by_cluster: Optional[Dict[str, Dict[str, MagnetosomeHit]]] = None,
    cwd: Optional[Path] = None,
    raise_on_error: bool = False,
) -> bool:
    """Run Clinker to generate comparison figure, then post-process colors and labels.

    Args:
        cwd: working directory for Clinker; relative output paths given in
            ``clinker_opts`` (e.g. ``-mo matrix.tsv``) are created there.
        raise_on_error: raise :class:`ClinkerError` instead of returning False.
    """
    if not gbk_files:
        return False

    def _fail(msg: str) -> bool:
        if raise_on_error:
            raise ClinkerError(msg)
        log.warning(msg)
        return False

    try:
        extra_args = parse_clinker_opts(clinker_opts)
    except ValueError as e:
        return _fail(str(e))

    clinker_bin = find_clinker()
    if not clinker_bin:
        return _fail("clinker is not installed (not found next to Python or on PATH); "
                     "interactive HTML figure was not generated.")

    plot_html = Path(plot_html).resolve()
    cmd = [
        clinker_bin,
        *[str(Path(p).resolve()) for p in gbk_files],
        "-p", str(plot_html),
        "-j", str(jobs),
        "-i", str(identity),
        *extra_args,
    ]

    if plot_html.exists():
        plot_html.unlink()
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True,
                       cwd=str(cwd) if cwd else None)
    except FileNotFoundError:
        return _fail("clinker is not installed on PATH; interactive HTML figure was not generated.")
    except subprocess.CalledProcessError as e:
        lines = [l.strip() for l in (e.stderr or e.stdout or "").splitlines() if l.strip()]
        err_lines = [l for l in lines if "error" in l.lower()] or lines[-1:]
        reason = f": {err_lines[-1]}" if err_lines else ""
        return _fail(f"clinker failed (exit code {e.returncode}){reason}. "
                     f"{plot_html.name} was not generated.")

    if not plot_html.exists():
        return _fail("clinker finished but did not write the HTML figure "
                     f"({plot_html.name}); check --clinker-opts.")
    # Apply custom coloring & protein labels
    customize_clinker_html(plot_html, hits, hits_by_cluster=hits_by_cluster)
    return True


def split_gbk_batches(
    gbk_files: List[Path],
    batch_size: int = DEFAULT_CLINKER_BATCH_SIZE,
) -> List[List[Path]]:
    """Partition *gbk_files* into sequential batches preserving input file order.

    ``batch_size <= 0`` disables splitting (all files are placed in a single batch).
    """
    files = list(gbk_files)
    if not files:
        return []
    if batch_size <= 0 or len(files) <= batch_size:
        return [files]
    return [files[i : i + batch_size] for i in range(0, len(files), batch_size)]


def batch_html_paths(base_html: Path, n_batches: int) -> List[Path]:
    """Return output HTML path(s) for *n_batches*.

    A single batch returns ``[base_html]`` (e.g. ``MGC_align.html``); multiple
    batches return zero-padded ``MGC_align_batch01.html``, ``MGC_align_batch02.html``, ...
    """
    if n_batches <= 0:
        return []
    if n_batches == 1:
        return [base_html]
    width = max(2, len(str(n_batches)))
    suffix = base_html.suffix or ".html"
    stem = base_html.stem if base_html.suffix else base_html.name
    parent = base_html.parent
    return [
        parent / f"{stem}_batch{idx:0{width}d}{suffix}"
        for idx in range(1, n_batches + 1)
    ]


def _with_batch_tag(raw_path: str, tag: str) -> str:
    p = Path(raw_path)
    if p.suffix:
        return str(p.with_name(f"{p.stem}_{tag}{p.suffix}"))
    return f"{raw_path}_{tag}"


def _per_batch_clinker_opts(
    clinker_opts: Optional[str],
    batch_idx: int,
    n_batches: int,
) -> Optional[str]:
    """Tag ``-o``/``--output`` and ``-mo``/``--matrix_out`` paths per batch so batches do not overwrite."""
    if not clinker_opts or n_batches <= 1:
        return clinker_opts
    tokens = parse_clinker_opts(clinker_opts)
    width = max(2, len(str(n_batches)))
    tag = f"batch{batch_idx:0{width}d}"
    file_flags = {"-o", "--output", "-mo", "--matrix_out"}
    out_tokens: List[str] = []
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok in file_flags and i + 1 < len(tokens):
            out_tokens.append(tok)
            out_tokens.append(_with_batch_tag(tokens[i + 1], tag))
            i += 2
        elif any(tok.startswith(f"{fl}=") for fl in ("--output", "--matrix_out")):
            fl, val = tok.split("=", 1)
            out_tokens.append(f"{fl}={_with_batch_tag(val, tag)}")
            i += 1
        else:
            out_tokens.append(tok)
            i += 1
    return shlex.join(out_tokens)


def run_clinker_batches(
    gbk_files: List[Path],
    base_html: Path,
    hits: Dict[str, MagnetosomeHit],
    batch_size: int = DEFAULT_CLINKER_BATCH_SIZE,
    jobs: int = 0,
    identity: float = 0.3,
    clinker_opts: Optional[str] = None,
    hits_by_cluster: Optional[Dict[str, Dict[str, MagnetosomeHit]]] = None,
    cwd: Optional[Path] = None,
    raise_on_error: bool = False,
    on_batch_start: Optional[Callable[[int, int, Path, int], None]] = None,
) -> List[Path]:
    """Run Clinker across one or more batches of GenBank files and return written HTML paths."""
    batches = split_gbk_batches(gbk_files, batch_size=batch_size)
    if not batches:
        return []
    html_targets = batch_html_paths(base_html, len(batches))
    written: List[Path] = []
    n_batches = len(batches)

    for idx, (batch_gbks, target_html) in enumerate(zip(batches, html_targets), 1):
        if on_batch_start is not None:
            on_batch_start(idx, n_batches, target_html, len(batch_gbks))
        batch_opts = _per_batch_clinker_opts(clinker_opts, idx, n_batches)
        ok = run_clinker_visualization(
            gbk_files=batch_gbks,
            plot_html=target_html,
            hits=hits,
            jobs=jobs,
            identity=identity,
            clinker_opts=batch_opts,
            hits_by_cluster=hits_by_cluster,
            cwd=cwd,
            raise_on_error=raise_on_error,
        )
        if ok:
            written.append(target_html)
    return written
