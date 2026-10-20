"""Clinker-compatible GenBank slice exporter and automated visualization runner."""

import json
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from ..core.cluster import MGCCluster
from ..core.predict import PredictedGene
from ..core.search import MagnetosomeHit
from ..utils.logger import log

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

        rec_id = f"{cluster.genome_name}_loc{i}" if len(slices) > 1 else cluster.cluster_id
        # GenBank LOCUS names/IDs must not contain whitespace
        rec_id = re.sub(r"\s+", "_", rec_id)
        rec_name = rec_id[:16]
        desc_label = f"slice {i}" if len(slices) > 1 else "cluster"
        new_rec = SeqRecord(
            sub_seq,
            id=rec_id,
            name=rec_name,
            description=f"Magnetosome Gene Cluster {desc_label} from {cluster.genome_name} ({ctg_id}:{offset+1}-{end_idx})"
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
