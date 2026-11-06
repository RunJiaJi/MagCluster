"""Tests for MagCluster v1.1.0 features and bug fixes."""

import gzip
import json
import re
import shutil
import sys
import warnings
from pathlib import Path

import pytest
from Bio import BiopythonWarning, SeqIO
from typer.testing import CliRunner

from magcluster import __version__
from magcluster.cli import app
from magcluster.core.cluster import MGCCluster, MGCSubcluster
from magcluster.core.predict import PredictedGene
from magcluster.core.search import MagnetosomeHit
from magcluster.viz.clinker_export import (
    DEFAULT_CLINKER_BATCH_SIZE,
    MAX_GBK_ID_LEN,
    MAX_LOCUS_LEN,
    _per_batch_clinker_opts,
    batch_html_paths,
    export_cluster_to_genbank,
    extract_hits_from_genbank,
    shorten_gbk_record_ids,
    split_gbk_batches,
)

DATA = Path(__file__).resolve().parent / "data"
MAI_GZ = DATA / "MSR-1_MAI90kb.fna.gz"
EXPECTED_GBK = DATA / "expected_1.0.0" / "MSR-1_MAI90kb_MGC.gbk"

runner = CliRunner()
HAS_CLINKER = shutil.which("clinker") is not None or (Path(sys.executable).parent / "clinker").exists()

ANSI_RE = re.compile(r"\x1b\[[0-9;?]*[ -/]*[@-~]")


def flat(text: str) -> str:
    return " ".join(ANSI_RE.sub("", text).split())


def run_cli(*args):
    res = runner.invoke(app, [str(a) for a in args])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert "Traceback" not in res.output
    return res


# --------------------------------------------------------------------------
# 1. Version bump
# --------------------------------------------------------------------------
def test_version_is_current():
    assert __version__ == "1.1.1"
    res = run_cli("--version")
    assert res.exit_code == 0 and "1.1.1" in flat(res.output)


# --------------------------------------------------------------------------
# 2. Long genome names: GenBank header keeps accession & emits zero warnings
# --------------------------------------------------------------------------
def test_long_genome_name_with_gca_accession_emits_no_biopython_warning(tmp_path):
    long_genome = "GCA_002753255.1_Candidatus_Lambdaproteobacteria_bacterium_PCRbin3"
    cluster_id = f"{long_genome}_MGC"

    g1 = PredictedGene("g1", "PDZZ01000224.1", 10, 909, 1, 0, "M" * 300)
    g2 = PredictedGene("g2", "PDZZ01000249.1", 3, 602, 1, 0, "M" * 200)
    g3 = PredictedGene("g3", "PDZZ01000185.1", 100, 7681, -1, 0, "M" * 200)

    sub1 = MGCSubcluster(f"{long_genome}_loc1", "PDZZ01000224.1", 10, 909, 10, 909, all_genes=[g1])
    sub2 = MGCSubcluster(f"{long_genome}_loc2", "PDZZ01000249.1", 3, 602, 3, 602, all_genes=[g2])
    sub3 = MGCSubcluster(f"{long_genome}_loc3", "PDZZ01000185.1", 100, 7681, 100, 7681, all_genes=[g3])

    cluster = MGCCluster(
        cluster_id=cluster_id,
        genome_name=long_genome,
        contig_id="PDZZ01000224.1,PDZZ01000249.1,PDZZ01000185.1",
        start_bp=3,
        end_bp=7681,
        cluster_core_start=3,
        cluster_core_end=7681,
        all_genes=[g1, g2, g3],
        subclusters=[sub1, sub2, sub3],
    )
    hits = {
        "g1": MagnetosomeHit("g1", "PDZZ01000224.1", "MamK", "MamK", 400.0, 1e-120, "core", "Other", "#EA4335", 10, 909, 1, 0, "M" * 300),
    }

    gbk_out = tmp_path / f"{cluster_id}.gbk"
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", BiopythonWarning)
        export_cluster_to_genbank(cluster, {}, hits, gbk_out)

    assert len(caught) == 0, [str(w.message) for w in caught]
    assert gbk_out.name == f"{cluster_id}.gbk"

    recs = list(SeqIO.parse(str(gbk_out), "genbank"))
    assert [r.id for r in recs] == [
        "GCA_002753255.1_loc1",
        "GCA_002753255.1_loc2",
        "GCA_002753255.1_loc3",
    ]
    # LOCUS names stay <= 16 chars and preserve _loc1/_loc2/_loc3
    assert all(len(r.name) <= MAX_LOCUS_LEN for r in recs)
    assert [r.name[-5:] for r in recs] == ["_loc1", "_loc2", "_loc3"]
    # Full genome name is still preserved in DEFINITION and /organism
    assert long_genome in recs[0].description
    src_feat = [f for f in recs[0].features if f.type == "source"][0]
    assert src_feat.qualifiers["organism"][0].replace(" ", "") == long_genome


def test_long_custom_genome_name_without_accession_never_warns(tmp_path):
    custom_long = "DeepSea_ColdSeep_Sample2026_Metagenome_Assembly_Bin_042_Candidatus_Magnetobacterium_sp_strain_X"
    rec_id, rec_name = shorten_gbk_record_ids(custom_long, f"{custom_long}_MGC", slice_idx=2, n_slices=3)
    assert len(rec_id) <= MAX_GBK_ID_LEN
    assert rec_id.endswith("_loc2")
    assert len(rec_name) <= MAX_LOCUS_LEN
    assert rec_name.endswith("_loc2")

    g1 = PredictedGene("g1", "c1", 1, 300, 1, 0, "M" * 100)
    cluster = MGCCluster(
        cluster_id=f"{custom_long}_MGC",
        genome_name=custom_long,
        contig_id="c1",
        start_bp=1,
        end_bp=300,
        cluster_core_start=1,
        cluster_core_end=300,
        all_genes=[g1],
    )
    gbk_out = tmp_path / f"{custom_long}_MGC.gbk"
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", BiopythonWarning)
        export_cluster_to_genbank(cluster, {}, {}, gbk_out)
    assert len(caught) == 0, [str(w.message) for w in caught]


def test_cli_run_long_filename_no_biopython_warning(tmp_path):
    long_fna = tmp_path / "GCA_002753255.1_Candidatus_Lambdaproteobacteria_bacterium_PCRbin3.fa"
    with gzip.open(MAI_GZ, "rt") as src:
        long_fna.write_text(src.read())
    out = tmp_path / "out"
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", BiopythonWarning)
        res = run_cli("run", long_fna, "-o", out, "--no-plot")
    assert res.exit_code == 0
    assert len(caught) == 0, [str(w.message) for w in caught]
    expected_gbk = out / "clusters_gbk" / "GCA_002753255.1_Candidatus_Lambdaproteobacteria_bacterium_PCRbin3_MGC.gbk"
    assert expected_gbk.exists()
    text = expected_gbk.read_text()
    assert "ACCESSION   GCA_002753255.1_loc1" in text


# --------------------------------------------------------------------------
# 3. Batching logic (default batch size = 50, input file order preserved)
# --------------------------------------------------------------------------
def test_default_batch_size_and_order_preserved():
    assert DEFAULT_CLINKER_BATCH_SIZE == 50
    paths = [Path(f"g_{i:03d}_MGC.gbk") for i in range(120)]
    batches = split_gbk_batches(paths)
    assert [len(b) for b in batches] == [50, 50, 20]
    assert [p for b in batches for p in b] == paths

    # <= batch_size stays a single batch named MGC_align.html
    assert len(split_gbk_batches(paths[:50])) == 1
    assert batch_html_paths(Path("/out/MGC_align.html"), 1) == [Path("/out/MGC_align.html")]
    assert batch_html_paths(Path("/out/MGC_align.html"), 3) == [
        Path("/out/MGC_align_batch01.html"),
        Path("/out/MGC_align_batch02.html"),
        Path("/out/MGC_align_batch03.html"),
    ]
    # batch_size=0 disables splitting
    assert len(split_gbk_batches(paths, batch_size=0)) == 1


def test_per_batch_clinker_opts_tags_output_files():
    assert _per_batch_clinker_opts("-na -mo matrix.tsv", 1, 1) == "-na -mo matrix.tsv"
    tagged = _per_batch_clinker_opts("-na -mo matrix.tsv -o align.txt", 2, 3)
    assert "matrix_batch02.tsv" in tagged
    assert "align_batch02.txt" in tagged


@pytest.mark.skipif(not HAS_CLINKER, reason="clinker not installed")
def test_run_and_clinker_batched_html_generation(tmp_path):
    in_dir = tmp_path / "genomes"
    in_dir.mkdir()
    for name in ("01_alpha.fna.gz", "02_beta.fna.gz", "03_gamma.fna.gz"):
        shutil.copy(MAI_GZ, in_dir / name)

    out_dir = tmp_path / "res"
    res = run_cli("run", in_dir, "-o", out_dir, "--batch-size", "2", "--clinker-opts", "-na -ufo")
    assert res.exit_code == 0
    b1 = out_dir / "MGC_align_batch01.html"
    b2 = out_dir / "MGC_align_batch02.html"
    assert b1.exists() and b2.exists()
    assert not (out_dir / "MGC_align.html").exists()

    # Check input order preserved across batches: 01 & 02 in batch01, 03 in batch02
    d1 = json.loads(re.search(r"<script>const data=(.*?);function serialise", b1.read_text()).group(1))
    d2 = json.loads(re.search(r"<script>const data=(.*?);function serialise", b2.read_text()).group(1))
    assert [c["name"] for c in d1["clusters"]] == ["01_alpha_MGC", "02_beta_MGC"]
    assert [c["name"] for c in d2["clusters"]] == ["03_gamma_MGC"]


# --------------------------------------------------------------------------
# 4. Standalone `magcluster clinker` command & help documentation
# --------------------------------------------------------------------------
def test_clinker_subcommand_help_and_options():
    top_help = flat(run_cli("--help").output)
    assert "clinker" in top_help and "run" in top_help and "db" in top_help

    clk_help = flat(run_cli("clinker", "--help").output)
    assert "Usage: magcluster clinker [OPTIONS]" in clk_help
    for opt in ("--batch-size", "--identity", "--list", "--clinker-opts", "--clinker-help", "--force"):
        assert opt in clk_help

    guide = flat(run_cli("clinker", "--clinker-help").output)
    assert "Recommended Clinker Options" in guide and "[default: 2]" in guide


def test_extract_hits_from_magcluster_gbk():
    hits = extract_hits_from_genbank(EXPECTED_GBK)
    assert len(hits) >= 25
    profiles = {h.profile for h in hits.values()}
    assert "MamK" in profiles and "MamE_full" in profiles


@pytest.mark.skipif(not HAS_CLINKER, reason="clinker not installed")
def test_standalone_clinker_command_workflows(tmp_path):
    # 1. Screen without plotting
    screen_out = tmp_path / "screen"
    assert run_cli("run", MAI_GZ, "-o", screen_out, "--no-plot").exit_code == 0
    assert not (screen_out / "MGC_align.html").exists()

    # 2. Pass the MagCluster results folder directly to `magcluster clinker` and write into it
    res1 = run_cli("clinker", screen_out, "-o", screen_out, "--clinker-opts", "-na")
    assert res1.exit_code == 0
    html_path = screen_out / "MGC_align.html"
    assert html_path.exists()

    # Verify MagCluster custom family colors and stripped _full labels were applied!
    m = re.search(r"<script>const data=(.*?);function serialise", html_path.read_text())
    data = json.loads(m.group(1))
    group_labels = {g["label"] for g in data["groups"]}
    assert "mam" in group_labels
    all_gene_labels = {
        g["label"]
        for c in data["clusters"]
        for loc in c["loci"]
        for g in loc["genes"]
        if g.get("is_mag")
    }
    assert "MamE" in all_gene_labels and "MamE_full" not in all_gene_labels

    # 3. Re-running without --force refuses to overwrite HTML, and --force keeps TSVs/GBKs intact!
    res_no_force = run_cli("clinker", screen_out, "-o", screen_out, "--clinker-opts", "-na")
    assert res_no_force.exit_code == 1 and "--force" in flat(res_no_force.output)

    res_force = run_cli("clinker", screen_out, "-o", screen_out, "--clinker-opts", "-na", "--force")
    assert res_force.exit_code == 0
    assert (screen_out / "MSR-1_MAI90kb_mgc_genes.tsv").exists()
    assert (screen_out / "mgc_summary.tsv").exists()
    assert (screen_out / "clusters_gbk" / "MSR-1_MAI90kb_MGC.gbk").exists()

    # 4. Direct .html target + --list ordering + .gbk.gz support
    gz_dir = tmp_path / "gz_gbks"
    gz_dir.mkdir()
    for tag in ("A_strain", "B_strain"):
        with open(EXPECTED_GBK, "rb") as src, gzip.open(gz_dir / f"{tag}_MGC.gbk.gz", "wb") as dst:
            dst.write(src.read())

    order_txt = tmp_path / "order.txt"
    order_txt.write_text("# Reverse order\nB_strain\nA_strain_MGC\n", encoding="utf-8")
    custom_html = tmp_path / "plots" / "custom_pair.html"
    res_list = run_cli("clinker", gz_dir, "-l", order_txt, "-o", custom_html, "--clinker-opts", "-na -ufo")
    assert res_list.exit_code == 0 and custom_html.exists()
    d_list = json.loads(re.search(r"<script>const data=(.*?);function serialise", custom_html.read_text()).group(1))
    assert [c["name"] for c in d_list["clusters"]] == ["B_strain_MGC", "A_strain_MGC"]


# --------------------------------------------------------------------------
# 5. Cohort summary table (mgc_summary.tsv) & multi-contig coordinate formatting
# --------------------------------------------------------------------------
def test_multi_contig_coordinate_formatting():
    g1 = PredictedGene("g1", "PDZZ01000224.1", 120, 5400, 1, 0, "M" * 100)
    g2 = PredictedGene("g2", "PDZZ01000249.1", 3, 4120, 1, 0, "M" * 100)
    g3 = PredictedGene("g3", "PDZZ01000185.1", 1100, 7681, 1, 0, "M" * 100)
    s1 = MGCSubcluster("x_loc1", "PDZZ01000224.1", 120, 5400, 120, 5400, all_genes=[g1])
    s2 = MGCSubcluster("x_loc2", "PDZZ01000249.1", 3, 4120, 3, 4120, all_genes=[g2])
    s3 = MGCSubcluster("x_loc3", "PDZZ01000185.1", 1100, 7681, 1100, 7681, all_genes=[g3])
    c = MGCCluster(
        cluster_id="x_MGC",
        genome_name="x",
        contig_id="PDZZ01000224.1,PDZZ01000249.1,PDZZ01000185.1",
        start_bp=3,
        end_bp=7681,
        cluster_core_start=3,
        cluster_core_end=7681,
        all_genes=[g1, g2, g3],
        subclusters=[s1, s2, s3],
    )
    assert c.format_coordinates() == (
        "PDZZ01000224.1:120-5400, PDZZ01000249.1:3-4120, PDZZ01000185.1:1100-7681"
    )
    assert c.format_coordinates(sep=";") == (
        "PDZZ01000224.1:120-5400;PDZZ01000249.1:3-4120;PDZZ01000185.1:1100-7681"
    )


def test_cohort_summary_table_written(tmp_path):
    neg_fna = tmp_path / "neg.fna"
    neg_fna.write_text(">ctg1\n" + ("ATG" + "GCT" * 60 + "TAA") + "\n", encoding="utf-8")
    out = tmp_path / "out"
    res = run_cli("run", MAI_GZ, neg_fna, "-o", out, "--no-plot")
    assert res.exit_code == 0

    summary = out / "mgc_summary.tsv"
    assert summary.exists()
    lines = summary.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 3
    header = lines[0].split("\t")
    row_pos = dict(zip(header, lines[1].split("\t")))
    row_neg = dict(zip(header, lines[2].split("\t")))

    assert row_pos["genome"] == "MSR-1_MAI90kb"
    assert row_pos["status"] == "MGC_detected"
    assert row_pos["mineral_type"] == "Magnetite"
    assert int(row_pos["n_distinct_families"]) >= 20
    assert row_pos["gbk_file"] == "clusters_gbk/MSR-1_MAI90kb_MGC.gbk"

    assert row_neg["genome"] == "neg"
    assert row_neg["status"] in ("no_MGC", "no_genes")
    assert row_neg["mineral_type"] == "-"


@pytest.mark.skipif(not HAS_CLINKER, reason="clinker not installed")
def test_external_genbank_without_magcluster_notes_is_annotated_and_colored(tmp_path):
    """External .gbk files lacking MagCluster /note metadata still get HMM-annotated & colored."""
    ext_gbk = tmp_path / "external_locus.gbk"
    recs = list(SeqIO.parse(str(EXPECTED_GBK), "genbank"))
    for rec in recs:
        rec.description = "External RefSeq genomic region"
        for idx_f, feat in enumerate(rec.features):
            if feat.type == "CDS":
                feat.qualifiers.pop("note", None)
                feat.qualifiers.pop("gene", None)
                feat.qualifiers.pop("label", None)
                feat.qualifiers.pop("protein", None)
                if idx_f % 2 == 0:
                    feat.qualifiers.pop("locus_tag", None)
                feat.qualifiers["product"] = ["hypothetical protein"]
    SeqIO.write(recs, str(ext_gbk), "genbank")

    out_html = tmp_path / "external_plot.html"
    res = run_cli("clinker", ext_gbk, "-o", out_html, "--clinker-opts", "-na")
    assert res.exit_code == 0 and out_html.exists()
    data = json.loads(re.search(r"<script>const data=(.*?);function serialise", out_html.read_text()).group(1))
    group_labels = {g["label"] for g in data["groups"]}
    assert "mam" in group_labels

