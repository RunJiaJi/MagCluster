"""Clinker-compatible GenBank slice exporter and automated visualization runner."""

import json
import re
import subprocess
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

def export_cluster_to_genbank(
    cluster: MGCCluster,
    contigs_dict: Dict[str, SeqRecord],
    hits: Dict[str, MagnetosomeHit],
    out_path: Path
):
    """Export sliced GenBank file for an MGC cluster (supporting single or multi-contig)."""
    # Group genes in this cluster by contig
    genes_by_ctg: Dict[str, List[PredictedGene]] = {}
    for g in cluster.all_genes:
        genes_by_ctg.setdefault(g.contig_id, []).append(g)

    records_to_write = []

    for ctg_id, ctg_genes in genes_by_ctg.items():
        ctg_genes.sort(key=lambda x: x.start)
        offset = ctg_genes[0].start - 1  # 0-based slice start
        end_idx = ctg_genes[-1].end      # 0-based slice end (exclusive)

        ctg_rec = contigs_dict.get(ctg_id)
        if ctg_rec and len(ctg_rec.seq) >= end_idx:
            sub_seq = ctg_rec.seq[offset:end_idx]
        else:
            span_len = max(1, end_idx - offset)
            sub_seq = Seq("N" * span_len)

        rec_name = cluster.cluster_id[:16]
        new_rec = SeqRecord(
            sub_seq,
            id=cluster.cluster_id,
            name=rec_name,
            description=f"Magnetosome Gene Cluster from {cluster.genome_name} ({ctg_id}:{offset+1}-{end_idx})"
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
                protein_name = h.profile
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
    out_dir: Path
) -> List[Path]:
    """Export sliced GenBank files for clinker visualization without creating extra CSV tables."""
    gbk_dir = out_dir / "clusters_gbk"
    gbk_dir.mkdir(parents=True, exist_ok=True)
    for old_gbk in gbk_dir.glob("*.gbk"):
        old_gbk.unlink()

    gbk_files = []
    for c in clusters:
        gbk_file = gbk_dir / f"{c.cluster_id}.gbk"
        genome_contigs = contigs_lookup.get(c.genome_name, {})
        export_cluster_to_genbank(c, genome_contigs, hits, gbk_file)
        gbk_files.append(gbk_file)

    return gbk_files

def customize_clinker_html(plot_html: Path, hits: Dict[str, MagnetosomeHit]) -> bool:
    """Post-process Clinker HTML to apply v0.2.12 family colors and protein name labels."""
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
        for locus in cluster.get("loci", []):
            for gene in locus.get("genes", []):
                names = gene.get("names", {})
                lt = names.get("locus_tag")
                gn = names.get("gene")

                hit = None
                if lt and lt in hits:
                    hit = hits[lt]
                elif gn and gn in hits:
                    hit = hits[gn]

                if hit:
                    protein_name = hit.profile
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

def run_clinker_visualization(
    gbk_files: List[Path],
    plot_html: Path,
    hits: Dict[str, MagnetosomeHit],
    jobs: int = 0,
    identity: float = 0.3,
    clinker_opts: Optional[str] = None
) -> bool:
    """Run Clinker to generate comparison figure, then post-process colors and labels."""
    if not gbk_files:
        return False

    import shlex

    cmd = [
        "clinker",
        *[str(p) for p in gbk_files],
        "-p", str(plot_html),
        "-j", str(jobs),
        "-i", str(identity)
    ]

    if clinker_opts:
        extra_args = shlex.split(clinker_opts)
        cmd.extend(extra_args)

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Apply custom coloring & protein labels
        customize_clinker_html(plot_html, hits)
        return True
    except FileNotFoundError:
        log.warning("clinker is not installed on PATH; skipped interactive HTML figure generation.")
        return False
    except subprocess.CalledProcessError as e:
        log.warning(f"clinker execution warning: {e.stderr}")
        return False
