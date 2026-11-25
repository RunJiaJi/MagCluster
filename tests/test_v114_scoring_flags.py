"""Tests for MagCluster v1.1.4 Contig_Edge, Multi_Contig, and 5-class mineral_type rules."""

import re
from pathlib import Path
from typing import Dict, List, Tuple

import pytest
from typer.testing import CliRunner

from magcluster.cli import app
from magcluster.core.cluster import ClusterEngine, MGCCluster, MGCSubcluster
from magcluster.core.predict import PredictedGene
from magcluster.core.scoring import ClusterScorer
from magcluster.core.search import MagnetosomeHit

runner = CliRunner()
ANSI_RE = re.compile(r"\x1b\[[0-9;?]*[ -/]*[@-~]")


def _make_hit(
    gid: str,
    ctg: str,
    profile: str,
    parent_family: str,
    role: str = "core",
    lineage: str = "Alpha_classical",
    start: int = 100,
    idx: int = 0,
) -> MagnetosomeHit:
    return MagnetosomeHit(
        gene_id=gid,
        contig_id=ctg,
        profile=profile,
        parent_family=parent_family,
        score=200.0,
        evalue=1e-40,
        role=role,
        lineage=lineage,
        color="#E63946",
        start=start,
        end=start + 899,
        strand=1,
        gene_index=idx,
        translation="M" * 200,
    )


def _build_contig(
    ctg: str,
    n_genes: int,
    hit_specs: Dict[int, Tuple[str, str, str, str]],
) -> Tuple[List[PredictedGene], Dict[str, MagnetosomeHit]]:
    """Build synthetic genes along one contig.

    ``hit_specs`` maps 0-based gene index -> (profile, parent_family, role, lineage).
    """
    genes: List[PredictedGene] = []
    hits: Dict[str, MagnetosomeHit] = {}
    for i in range(n_genes):
        gid = f"{ctg}_g{i}"
        st = 1 + i * 1000
        genes.append(PredictedGene(gid, ctg, st, st + 899, 1, i, "M" * 200))
        if i in hit_specs:
            prof, fam, role, lineage = hit_specs[i]
            hits[gid] = _make_hit(gid, ctg, prof, fam, role=role, lineage=lineage, start=st, idx=i)
    return genes, hits


def _make_cluster_from_hits(hits: List[MagnetosomeHit], genome_name: str = "sp") -> MGCCluster:
    core = [h for h in hits if h.role == "core"]
    acc = [h for h in hits if h.role == "accessory"]
    iron = [h for h in hits if h.role == "iron_context"]
    genes = [
        PredictedGene(h.gene_id, h.contig_id, h.start, h.end, h.strand, i, h.translation)
        for i, h in enumerate(hits)
    ]
    sub = MGCSubcluster(
        subcluster_id=f"{genome_name}_loc1",
        contig_id="ctg1",
        start_bp=1,
        end_bp=10000,
        cluster_core_start=1,
        cluster_core_end=10000,
        core_hits=core,
        accessory_hits=acc,
        iron_hits=iron,
        all_genes=genes,
        on_edge=False,
    )
    return MGCCluster(
        cluster_id=f"{genome_name}_MGC",
        genome_name=genome_name,
        contig_id="ctg1",
        start_bp=1,
        end_bp=10000,
        cluster_core_start=1,
        cluster_core_end=10000,
        core_hits=core,
        accessory_hits=acc,
        iron_hits=iron,
        all_genes=genes,
        subclusters=[sub],
        on_edge=False,
    )


# --------------------------------------------------------------------------
# 1. Contig_Edge decoupled from --flank & circular replicons excluded
# --------------------------------------------------------------------------
@pytest.mark.parametrize("flank", [0, 3, 10])
def test_contig_edge_decoupled_from_flank_interior_never_flagged(flank):
    """Magnetosome proteins at index 1..N-2 on a linear contig are never Contig_Edge regardless of --flank."""
    genes, hits = _build_contig(
        "ctg1",
        10,
        {
            1: ("MamA", "MamA", "core", "Alpha_classical"),
            2: ("MamB", "MamB", "core", "Alpha_classical"),
            3: ("MamE_full", "MamE", "core", "Alpha_classical"),
            4: ("MamK", "MamK", "core", "Alpha_classical"),
        },
    )
    clusters = ClusterEngine(flank_genes=flank).cluster_genome("sp", genes, hits)
    assert len(clusters) == 1
    c = clusters[0]
    assert c.subclusters[0].on_edge is False
    assert c.on_edge is False
    sc = ClusterScorer.score_cluster(c)
    assert sc.flags == []


@pytest.mark.parametrize("flank", [0, 3, 10])
@pytest.mark.parametrize("positions", [(0, 1, 2), (7, 8, 9)])
def test_contig_edge_flagged_when_magnetosome_gene_at_linear_boundary(flank, positions):
    """Magnetosome protein directly at index 0 or N-1 on a linear contig triggers Contig_Edge regardless of --flank."""
    p0, p1, p2 = positions
    genes, hits = _build_contig(
        "ctg1",
        10,
        {
            p0: ("MamA", "MamA", "core", "Alpha_classical"),
            p1: ("MamB", "MamB", "core", "Alpha_classical"),
            p2: ("MamE_full", "MamE", "core", "Alpha_classical"),
        },
    )
    clusters = ClusterEngine(flank_genes=flank).cluster_genome("sp", genes, hits)
    assert len(clusters) == 1
    c = clusters[0]
    assert c.subclusters[0].on_edge is True
    assert c.on_edge is True
    sc = ClusterScorer.score_cluster(c)
    assert sc.flags == ["Contig_Edge"]


def test_circular_replicon_never_flagged_as_contig_edge():
    """Complete circular contigs never get Contig_Edge even when genes sit at index 0 / N-1 or span origin."""
    # 1. Origin-spanning two-subcluster MGC on a circular chromosome (enough interior genes so flank=3 stays 2 slices)
    genes, hits = [], {}
    starts = (
        [1000 + i * 1000 for i in range(7)]           # g0..g6 near origin start (g0..g2 are hits)
        + [40_000 + i * 2000 for i in range(4)]       # g7..g10 mid-chromosome filler
        + [90_000 + i * 1000 for i in range(7)]       # g11..g17 near origin end (g15..g17 are hits)
    )
    hit_fams = {0: "MamA", 1: "MamB", 2: "MamE", 15: "MamK", 16: "MamM", 17: "MamO"}
    for i, st in enumerate(starts):
        genes.append(PredictedGene(f"g{i}", "chr", st, st + 899, 1, i, "M" * 200))
        if i in hit_fams:
            fam = hit_fams[i]
            hits[f"g{i}"] = _make_hit(f"g{i}", "chr", fam, fam, start=st, idx=i)

    engine = ClusterEngine(flank_genes=3)
    circ = engine.cluster_genome("sp", genes, hits, {"chr": 100_000}, {"chr": True})[0]
    assert len(circ.subclusters) == 2
    assert circ.wraps_origin is True
    assert all(not sub.on_edge for sub in circ.subclusters)
    assert circ.on_edge is False
    sc_circ = ClusterScorer.score_cluster(circ)
    assert sc_circ.flags == ["Spans_Origin"]

    # Same coordinates on a linear contig must report per-subcluster Contig_Edge(loc1;loc2)
    linear = engine.cluster_genome("sp", genes, hits, {"chr": 100_000}, {"chr": False})[0]
    assert linear.wraps_origin is False
    assert [sub.on_edge for sub in linear.subclusters] == [True, True]
    assert linear.on_edge is False  # no genome-wide any() aggregation for multi-subcluster
    sc_lin = ClusterScorer.score_cluster(linear)
    assert sc_lin.flags == ["Contig_Edge(loc1;loc2)"]

    # 2. Single-subcluster MGC starting at gene index 0 on a circular replicon -> still no Contig_Edge
    g_single, h_single = _build_contig(
        "chr_circ",
        10,
        {
            0: ("MamA", "MamA", "core", "Alpha_classical"),
            1: ("MamB", "MamB", "core", "Alpha_classical"),
            2: ("MamE_full", "MamE", "core", "Alpha_classical"),
        },
    )
    c_single = engine.cluster_genome("sp_circ", g_single, h_single, {"chr_circ": 20_000}, {"chr_circ": True})[0]
    assert c_single.on_edge is False
    assert ClusterScorer.score_cluster(c_single).flags == []


# --------------------------------------------------------------------------
# 2. Per-subcluster Contig_Edge(loc...) & Multi_Contig (no genome-wide any())
# --------------------------------------------------------------------------
def test_multi_subcluster_same_contig_reports_specific_edge_locus_only():
    """Two distant subclusters on one linear contig: only loc2 sits at the end -> Contig_Edge(loc2), no Multi_Contig."""
    genes, hits = _build_contig(
        "chr1",
        30,
        {
            5: ("MamA", "MamA", "core", "Alpha_classical"),
            6: ("MamB", "MamB", "core", "Alpha_classical"),
            7: ("MamE_full", "MamE", "core", "Alpha_classical"),
            27: ("MamK", "MamK", "core", "Alpha_classical"),
            28: ("MamM", "MamM", "core", "Alpha_classical"),
            29: ("MamO_full", "MamO", "core", "Alpha_classical"),
        },
    )
    c = ClusterEngine(flank_genes=3).cluster_genome("amb1", genes, hits)[0]
    assert len(c.subclusters) == 2
    assert [s.on_edge for s in c.subclusters] == [False, True]
    assert c.on_edge is False
    assert c.is_multi_contig is False
    sc = ClusterScorer.score_cluster(c)
    assert sc.flags == ["Contig_Edge(loc2)"]


def test_multi_contig_flags_and_specific_edge_subclusters():
    """Multi-contig MGC gets Multi_Contig and lists exact edge loci (loc1;loc3)."""
    g1, h1 = _build_contig("ctgA", 8, {
        0: ("MamA", "MamA", "core", "Alpha_classical"),
        1: ("MamB", "MamB", "core", "Alpha_classical"),
        2: ("MamE_full", "MamE", "core", "Alpha_classical"),
    })
    g2, h2 = _build_contig("ctgB", 10, {
        3: ("MamI", "MamI", "core", "Alpha_classical"),
        4: ("MamK", "MamK", "core", "Alpha_classical"),
        5: ("MamL", "MamL", "core", "Alpha_classical"),
    })
    g3, h3 = _build_contig("ctgC", 8, {
        5: ("MamM", "MamM", "core", "Alpha_classical"),
        6: ("MamO_full", "MamO", "core", "Alpha_classical"),
        7: ("MamQ", "MamQ", "core", "Alpha_classical"),
    })
    c = ClusterEngine(flank_genes=3).cluster_genome("mag1", g1 + g2 + g3, {**h1, **h2, **h3})[0]
    assert len(c.subclusters) == 3
    assert [s.on_edge for s in c.subclusters] == [True, False, True]
    assert c.on_edge is False
    assert c.is_multi_contig is True
    sc = ClusterScorer.score_cluster(c)
    assert sc.flags == ["Contig_Edge(loc1;loc3)", "Multi_Contig"]


def test_multi_contig_without_edge_reports_multi_contig_only():
    g1, h1 = _build_contig("ctgA", 8, {
        2: ("MamA", "MamA", "core", "Alpha_classical"),
        3: ("MamB", "MamB", "core", "Alpha_classical"),
        4: ("MamE_full", "MamE", "core", "Alpha_classical"),
    })
    g2, h2 = _build_contig("ctgB", 8, {
        2: ("MamI", "MamI", "core", "Alpha_classical"),
        3: ("MamK", "MamK", "core", "Alpha_classical"),
        4: ("MamM", "MamM", "core", "Alpha_classical"),
    })
    c = ClusterEngine(flank_genes=3).cluster_genome("mag2", g1 + g2, {**h1, **h2})[0]
    assert [s.on_edge for s in c.subclusters] == [False, False]
    assert c.is_multi_contig is True
    sc = ClusterScorer.score_cluster(c)
    assert sc.flags == ["Multi_Contig"]


# --------------------------------------------------------------------------
# 3. 5-class mineral_type classification
# --------------------------------------------------------------------------
def test_mineral_type_magnetite_and_unknown():
    # 3 distinct non-greigite core families (including MamK) -> Magnetite
    c_mag = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamA", "MamA", "core", "Alpha_classical"),
        _make_hit("g2", "ctg1", "MamB", "MamB", "core", "Alpha_classical"),
        _make_hit("g3", "ctg1", "MamK", "MamK", "core", "Alpha_classical"),
    ])
    sc_mag = ClusterScorer.score_cluster(c_mag)
    assert sc_mag.mineral_type == "Magnetite"
    assert sc_mag.greigite_families == []
    assert sc_mag.magnetite_core_families == ["MamA", "MamB"]

    # Accessory-only cluster (e.g. Mad1 + Mad2 + Mad10) -> Unknown
    c_acc = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "Mad1", "Mad1", "accessory", "Desulfo_Mad"),
        _make_hit("g2", "ctg1", "Mad2", "Mad2", "accessory", "Desulfo_Mad"),
        _make_hit("g3", "ctg1", "Mad10", "Mad10", "accessory", "Desulfo_Mad"),
    ])
    assert ClusterScorer.score_cluster(c_acc).mineral_type == "Unknown"

    # Only 2 non-greigite core families + accessories -> Unknown
    c_partial = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamB", "MamB", "core", "Alpha_classical"),
        _make_hit("g2", "ctg1", "MamK", "MamK", "core", "Alpha_classical"),
        _make_hit("g3", "ctg1", "Mad10", "Mad10", "accessory", "Desulfo_Mad"),
        _make_hit("g4", "ctg1", "Mad23", "Mad23", "accessory", "Desulfo_Mad"),
    ])
    assert ClusterScorer.score_cluster(c_partial).mineral_type == "Unknown"


def test_mineral_type_putative_greigite_and_family_deduplication():
    # Split MamE (MamE_Nter_greigite + MamE_Cter_greigite) deduplicates to G=1 -> Putative_Greigite
    c_split_mame = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamE_Nter_greigite", "MamE", "core", "Desulfo_greigite"),
        _make_hit("g2", "ctg1", "MamE_Cter_greigite", "MamE", "core", "Desulfo_greigite"),
        _make_hit("g3", "ctg1", "Mad1", "Mad1", "accessory", "Desulfo_Mad"),
        _make_hit("g4", "ctg1", "Mad2", "Mad2", "accessory", "Desulfo_Mad"),
    ])
    sc_split = ClusterScorer.score_cluster(c_split_mame)
    assert sc_split.mineral_type == "Putative_Greigite"
    assert sc_split.greigite_families == ["MamE"]

    # Tandem duplicated MamT_greigite (3 copies) still deduplicates to G=1 -> Putative_Greigite
    c_tandem = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamT_greigite", "MamT", "accessory", "Desulfo_greigite"),
        _make_hit("g2", "ctg1", "MamT_greigite", "MamT", "accessory", "Desulfo_greigite"),
        _make_hit("g3", "ctg1", "MamT_greigite", "MamT", "accessory", "Desulfo_greigite"),
        _make_hit("g4", "ctg1", "Mad1", "Mad1", "accessory", "Desulfo_Mad"),
        _make_hit("g5", "ctg1", "Mad2", "Mad2", "accessory", "Desulfo_Mad"),
    ])
    sc_tandem = ClusterScorer.score_cluster(c_tandem)
    assert sc_tandem.mineral_type == "Putative_Greigite"
    assert sc_tandem.greigite_families == ["MamT"]

    # 2 distinct _greigite families -> Putative_Greigite
    c_two = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamB_greigite", "MamB", "core", "Desulfo_greigite"),
        _make_hit("g2", "ctg1", "MamQ_greigite", "MamQ", "core", "Desulfo_greigite"),
        _make_hit("g3", "ctg1", "Mad10", "Mad10", "accessory", "Desulfo_Mad"),
    ])
    sc_two = ClusterScorer.score_cluster(c_two)
    assert sc_two.mineral_type == "Putative_Greigite"
    assert sc_two.greigite_families == ["MamB", "MamQ"]


def test_mineral_type_greigite_excludes_shared_mamk_from_magnetite_core():
    # Pure greigite MGC: 3 greigite families (incl. MamT_greigite) + shared MamK + 2 regular core (no dual copy)
    # Excluding MamK leaves M=2 (<3) and dual_copy=0 -> Greigite (not falsely Both!)
    c_greig = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamA_greigite", "MamA", "core", "Desulfo_greigite"),
        _make_hit("g2", "ctg1", "MamB_greigite", "MamB", "core", "Desulfo_greigite"),
        _make_hit("g3", "ctg1", "MamT_greigite", "MamT", "accessory", "Desulfo_greigite"),
        _make_hit("g4", "ctg1", "MamK", "MamK", "core", "Alpha_classical"),
        _make_hit("g5", "ctg1", "MamO_domain", "MamO", "core", "Alpha_classical"),
        _make_hit("g6", "ctg1", "MamP", "MamP", "core", "Alpha_classical"),
    ])
    sc = ClusterScorer.score_cluster(c_greig)
    assert sc.mineral_type == "Greigite"
    assert sc.greigite_families == ["MamA", "MamB", "MamT"]
    assert sc.magnetite_core_families == ["MamO", "MamP"]
    assert sc.dual_copy_families == []

    # Even with 1 dual copy, M=1 (<2) still stays Greigite
    c_greig_1m = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamA_greigite", "MamA", "core", "Desulfo_greigite"),
        _make_hit("g2", "ctg1", "MamB_greigite", "MamB", "core", "Desulfo_greigite"),
        _make_hit("g3", "ctg1", "MamQ_greigite", "MamQ", "core", "Desulfo_greigite"),
        _make_hit("g4", "ctg1", "MamB", "MamB", "core", "Alpha_classical"),
        _make_hit("g5", "ctg1", "MamK", "MamK", "core", "Alpha_classical"),
    ])
    sc_1m = ClusterScorer.score_cluster(c_greig_1m)
    assert sc_1m.mineral_type == "Greigite"
    assert sc_1m.dual_copy_families == ["MamB"]


def test_mineral_type_both_full_dual_and_dual_copy_cooccurrence():
    # Case A (BW-1 style): G >= 3 and M >= 3 (excluding MamK) -> Both
    c_bw1 = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamA_greigite", "MamA", "core", "Desulfo_greigite"),
        _make_hit("g2", "ctg1", "MamB_greigite", "MamB", "core", "Desulfo_greigite"),
        _make_hit("g3", "ctg1", "MamE_greigite", "MamE", "core", "Desulfo_greigite"),
        _make_hit("g4", "ctg1", "MamQ_greigite", "MamQ", "core", "Desulfo_greigite"),
        _make_hit("g5", "ctg1", "MamM", "MamM", "core", "Alpha_classical"),
        _make_hit("g6", "ctg1", "MamL", "MamL", "core", "Alpha_classical"),
        _make_hit("g7", "ctg1", "MamB", "MamB", "core", "Alpha_classical"),
        _make_hit("g8", "ctg1", "MamK", "MamK", "core", "Alpha_classical"),
    ])
    sc_bw1 = ClusterScorer.score_cluster(c_bw1)
    assert sc_bw1.mineral_type == "Both"
    assert sc_bw1.greigite_families == ["MamA", "MamB", "MamE", "MamQ"]
    assert sc_bw1.magnetite_core_families == ["MamB", "MamL", "MamM"]
    assert sc_bw1.dual_copy_families == ["MamB"]

    # Case B (fragmented dual-biomineralizing MAG): G >= 3, M == 2, AND dual-copy co-occurrence -> Both
    c_dual_copy = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamA_greigite", "MamA", "core", "Desulfo_greigite"),
        _make_hit("g2", "ctg1", "MamB_greigite", "MamB", "core", "Desulfo_greigite"),
        _make_hit("g3", "ctg1", "MamI_greigite", "MamI", "core", "Desulfo_greigite"),
        _make_hit("g4", "ctg1", "MamB", "MamB", "core", "Alpha_classical"),
        _make_hit("g5", "ctg1", "MamM", "MamM", "core", "Alpha_classical"),
    ])
    sc_dc = ClusterScorer.score_cluster(c_dual_copy)
    assert sc_dc.mineral_type == "Both"
    assert sc_dc.dual_copy_families == ["MamB"]
    assert sc_dc.magnetite_core_families == ["MamB", "MamM"]


def test_greigite_detection_uses_structured_lineage_metadata():
    """Profiles tagged lineage='Desulfo_greigite' are recognized even without 'greigite' in profile name."""
    c = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "CustomA", "MamA", "core", "Desulfo_greigite"),
        _make_hit("g2", "ctg1", "CustomB", "MamB", "core", "Desulfo_greigite"),
        _make_hit("g3", "ctg1", "CustomQ", "MamQ", "core", "Desulfo_greigite"),
    ])
    assert ClusterScorer.score_cluster(c).mineral_type == "Greigite"


def test_cli_summary_writes_new_flags_for_multi_contig_edge_genome(tmp_path):
    """End-to-end CLI test on a 2-contig Prodigal .faa verifying Multi_Contig & Contig_Edge(loc1) in mgc_summary.tsv."""
    from Bio import SeqIO

    gbk = Path(__file__).resolve().parent / "data" / "expected_1.0.0" / "MSR-1_MAI90kb_MGC.gbk"
    prots: Dict[str, str] = {}
    for rec in SeqIO.parse(str(gbk), "genbank"):
        for feat in rec.features:
            if feat.type == "CDS" and "gene" in feat.qualifiers and "translation" in feat.qualifiers:
                prots.setdefault(feat.qualifiers["gene"][0], feat.qualifiers["translation"][0])

    # Build two contigs:
    # ctg1: MamA at index 0 (edge!), MamB, MamE, non-mag filler
    # ctg2: non-mag filler, MamK, MamM, MamO, non-mag filler (interior -> not edge)
    filler = "M" + "ACDEFGHIKLMNPQRSTVWY" * 6
    lines = [
        f">ctg1_1 # 1 # 900 # 1 # ID=1_1;contig=ctg1\n{prots['MamA']}",
        f">ctg1_2 # 1001 # 1900 # 1 # ID=1_2;contig=ctg1\n{prots['MamB']}",
        f">ctg1_3 # 2001 # 3900 # 1 # ID=1_3;contig=ctg1\n{prots['MamE']}",
        f">ctg1_4 # 4001 # 4400 # 1 # ID=1_4;contig=ctg1\n{filler}",
        f">ctg2_1 # 1 # 400 # 1 # ID=2_1;contig=ctg2\n{filler}",
        f">ctg2_2 # 1001 # 2100 # 1 # ID=2_2;contig=ctg2\n{prots['MamK']}",
        f">ctg2_3 # 2201 # 3200 # 1 # ID=2_3;contig=ctg2\n{prots['MamM']}",
        f">ctg2_4 # 3301 # 5200 # 1 # ID=2_4;contig=ctg2\n{prots['MamO']}",
        f">ctg2_5 # 5301 # 5700 # 1 # ID=2_5;contig=ctg2\n{filler}",
    ]
    faa = tmp_path / "split_mag.faa"
    faa.write_text("\n".join(lines) + "\n", encoding="utf-8")

    out = tmp_path / "out"
    res = runner.invoke(app, ["run", str(faa), "-o", str(out), "--no-plot"])
    assert res.exit_code == 0, res.output

    summary = (out / "mgc_summary.tsv").read_text(encoding="utf-8").splitlines()
    row = dict(zip(summary[0].split("\t"), summary[1].split("\t")))
    assert row["status"] == "MGC_detected"
    assert row["mineral_type"] == "Magnetite"
    assert row["flags"] == "Contig_Edge(loc1),Multi_Contig"


# --------------------------------------------------------------------------
# 4. P0/P1 scientific correctness & edge-case regressions
# --------------------------------------------------------------------------
def test_flank_never_absorbs_unclustered_hits_into_mgc(tmp_path):
    """Increasing --flank must not pull unclustered magnetosome hits into MGC membership or mineral scoring."""
    from magcluster.utils.io import export_genome_magnetosome_fasta, write_genome_cluster_tables

    # Core cluster at indices 10, 11, 12 (MamA, MamB, MamE); distant unclustered MamQ_greigite at index 3
    # (gene_gap = 6 > max_gap = 5), which falls inside flank_genes=10!
    genes, hits = _build_contig(
        "ctg1",
        20,
        {
            3: ("MamQ_greigite", "MamQ", "core", "Desulfo_greigite"),
            10: ("MamA", "MamA", "core", "Alpha_classical"),
            11: ("MamB", "MamB", "core", "Alpha_classical"),
            12: ("MamE_full", "MamE", "core", "Alpha_classical"),
        },
    )
    for g in genes:
        g.genome_name = "sp"

    c = ClusterEngine(max_intergenic_genes=5, flank_genes=10).cluster_genome("sp", genes, hits)[0]
    # Gene 3 is in flank region_genes for visual context, but NOT in clustered hits or mineral scoring!
    assert "ctg1_g3" in [g.gene_id for g in c.all_genes]
    assert "ctg1_g3" not in c.clustered_hit_ids
    assert c.distinct_magnetosome_families == {"MamA", "MamB", "MamE"}
    assert ClusterScorer.score_cluster(c).mineral_type == "Magnetite"

    # Per-genome TSV and FAA must report ctg1_g3 as unclustered, not MGC
    lookup = {g.gene_id: g for g in genes}
    tsv = write_genome_cluster_tables("sp", [c], hits, tmp_path, genes_lookup=lookup)
    faa = export_genome_magnetosome_fasta("sp", [c], hits, tmp_path, genes_lookup=lookup)
    tsv_by_id = {r.split("\t")[1]: r.split("\t")[0] for r in tsv.read_text().splitlines()[1:]}
    assert tsv_by_id == {
        "ctg1_g10": "MGC",
        "ctg1_g11": "MGC",
        "ctg1_g12": "MGC",
        "ctg1_g3": "unclustered",
    }
    assert ">ctg1_g3 cluster=unclustered " in faa.read_text()


def test_native_prokka_locus_tags_matching_genome_name_are_not_corrupted(tmp_path):
    """GenBank locus_tags like AMB1_0001 must remain AMB1_0001 in TSV/FAA (matching GBK), not gene_0001."""
    from magcluster.utils.io import export_genome_magnetosome_fasta, write_genome_cluster_tables

    g1 = PredictedGene("AMB1_0001", "ctg1", 100, 999, 1, 0, "M" * 200, locus_tag="AMB1_0001")
    g1.genome_name = "AMB1"
    h1 = _make_hit("AMB1_0001", "ctg1", "MamA", "MamA", start=100, idx=0)
    c = _make_cluster_from_hits([h1], genome_name="AMB1")
    c.all_genes = [g1]

    tsv = write_genome_cluster_tables("AMB1", [c], {"AMB1_0001": h1}, tmp_path, genes_lookup={"AMB1_0001": g1})
    faa = export_genome_magnetosome_fasta("AMB1", [c], {"AMB1_0001": h1}, tmp_path, genes_lookup={"AMB1_0001": g1})
    assert tsv.read_text().splitlines()[1].split("\t")[1] == "AMB1_0001"
    assert faa.read_text().splitlines()[0].startswith(">AMB1_0001 cluster=MGC ")


def test_genbank_none_strand_and_origin_wrapping_compound_cds(tmp_path):
    """parse_genbank tolerates strand=None and anchors circular origin-wrapping CompoundLocation CDSs."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    from magcluster.core.parser import parse_genbank

    rec = SeqRecord(Seq("ATG" + "GCT" * 3332 + "TAA"), id="NC_CIRC.1", name="NC_CIRC")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular"

    # CDS 1: strand=None
    f_none = SeqFeature(
        FeatureLocation(1000, 1300, strand=None),
        type="CDS",
        qualifiers={"locus_tag": ["LT_NONE"], "translation": ["M" * 100]},
    )
    # CDS 2: origin-wrapping CompoundLocation join(9600..10002, 0..200) -> 602 bp total, not 10,002 bp!
    wrap_loc = CompoundLocation([FeatureLocation(9600, len(rec.seq), strand=1), FeatureLocation(0, 200, strand=1)])
    f_wrap = SeqFeature(
        wrap_loc,
        type="CDS",
        qualifiers={"locus_tag": ["LT_WRAP"], "translation": ["M" * 200]},
    )
    rec.features = [f_none, f_wrap]

    gbk = tmp_path / "circ.gbff"
    SeqIO.write([rec], str(gbk), "genbank")
    parsed = parse_genbank(gbk)
    by_id = {g.gene_id: g for g in parsed.genes}
    assert by_id["LT_NONE"].strand == 1
    assert by_id["LT_WRAP"].start == 9601
    assert by_id["LT_WRAP"].end == len(rec.seq)


def test_overlapping_nested_cds_not_truncated_in_cluster_or_gbk(tmp_path):
    """When a shorter CDS is nested inside a longer CDS, end_bp and sliced GBK preserve the outer CDS end."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from magcluster.viz.clinker_export import export_cluster_to_genbank

    genes = [
        PredictedGene("g0", "ctg1", 100, 900, 1, 0, "M" * 200),
        PredictedGene("g1", "ctg1", 1000, 1900, 1, 1, "M" * 200),
        PredictedGene("g2_outer", "ctg1", 2000, 4000, 1, 2, "M" * 500),
        PredictedGene("g3_nested", "ctg1", 2500, 2800, -1, 3, "M" * 90),
    ]
    hits = {
        "g0": _make_hit("g0", "ctg1", "MamA", "MamA", start=100, idx=0),
        "g1": _make_hit("g1", "ctg1", "MamB", "MamB", start=1000, idx=1),
        "g2_outer": _make_hit("g2_outer", "ctg1", "MamE_full", "MamE", start=2000, idx=2),
    }
    hits["g2_outer"].end = 4000

    c = ClusterEngine(flank_genes=1).cluster_genome("sp", genes, hits)[0]
    assert c.subclusters[0].end_bp == 4000
    assert c.subclusters[0].cluster_core_end == 4000

    ctg_rec = SeqRecord(Seq("A" * 5000), id="ctg1")
    gbk = tmp_path / "nested.gbk"
    export_cluster_to_genbank(c, {"ctg1": ctg_rec}, hits, gbk)
    out_rec = SeqIO.read(str(gbk), "genbank")
    assert len(out_rec.seq) == 4000 - 100 + 1


def test_mixed_faa_skips_unplaced_coordless_genes():
    """Mixed .faa with Prodigal headers + plain headers must never cluster coord-less unplaced hits."""
    genes = [
        PredictedGene("p1", "ctg1", 100, 900, 1, 0, "M" * 200),
        PredictedGene("u1", "unplaced", 0, 0, 0, 0, "M" * 200),
        PredictedGene("u2", "unplaced", 0, 0, 0, 1, "M" * 200),
        PredictedGene("u3", "unplaced", 0, 0, 0, 2, "M" * 200),
    ]
    hits = {
        "u1": _make_hit("u1", "unplaced", "MamA", "MamA", start=0, idx=0),
        "u2": _make_hit("u2", "unplaced", "MamB", "MamB", start=0, idx=1),
        "u3": _make_hit("u3", "unplaced", "MamE_full", "MamE", start=0, idx=2),
    }
    assert ClusterEngine().cluster_genome("mixed", genes, hits) == []


def test_linear_or_partial_complete_sequence_header_not_circular():
    from magcluster.core.parser import _is_circular_fasta_header

    assert _is_circular_fasta_header("NC_007626.1 Magnetospirillum magneticum AMB-1, complete genome") is True
    assert _is_circular_fasta_header("ctg1 topology=circular") is True
    assert _is_circular_fasta_header("AF123456 linear plasmid, complete sequence") is False
    assert _is_circular_fasta_header("ctg1 complete genome topology=linear") is False
    assert _is_circular_fasta_header("ctg1 complete sequence, partial") is False


def test_circular_origin_rescue_and_distant_islands_not_collapsed():
    """Circular replicon rescues a 2+2 MGC split across origin, while keeping 400 kb distant islands separate."""
    engine = ClusterEngine(max_intergenic_genes=5, max_distance_bp=25_000, min_distinct_families=3, flank_genes=1)

    # 1. Split MGC: 2 families at end + 2 families at start (each side < 3 alone, combined = 4 >= 3)
    starts = [500, 1600, 3000, 50_000, 96_000, 97_500, 98_600]
    fams = {0: "MamA", 1: "MamB", 5: "MamE", 6: "MamQ"}
    genes, hits = [], {}
    for i, st in enumerate(starts):
        genes.append(PredictedGene(f"g{i}", "chr", st, st + 899, 1, i, "M" * 200))
        if i in fams:
            hits[f"g{i}"] = _make_hit(f"g{i}", "chr", fams[i], fams[i], start=st, idx=i)

    # On a linear contig, neither side reaches 3 families -> no MGC
    assert engine.cluster_genome("sp", genes, hits, {"chr": 100_000}, {"chr": False}) == []
    # On a circular replicon, both halves bridge across 0 bp -> rescued & Spans_Origin!
    circ = engine.cluster_genome("sp", genes, hits, {"chr": 100_000}, {"chr": True})
    assert len(circ) == 1
    assert circ[0].wraps_origin is True
    assert circ[0].distinct_magnetosome_families == {"MamA", "MamB", "MamE", "MamQ"}
    assert ClusterScorer.score_cluster(circ[0]).flags == ["Spans_Origin"]

    # 2. Two distant islands on a 4 Mb circular chromosome (200 kb & 3800 kb -> 400 kb across origin > 25 kb)
    d_starts = [200_000, 201_000, 202_000, 3_800_000, 3_801_000, 3_802_000]
    d_fams = ["MamA", "MamB", "MamE", "MamK", "MamM", "MamO"]
    d_genes, d_hits = [], {}
    for i, (st, fam) in enumerate(zip(d_starts, d_fams)):
        d_genes.append(PredictedGene(f"d{i}", "chr", st, st + 899, 1, i, "M" * 200))
        d_hits[f"d{i}"] = _make_hit(f"d{i}", "chr", fam, fam, start=st, idx=i)

    c_dist = engine.cluster_genome("dist", d_genes, d_hits, {"chr": 4_000_000}, {"chr": True})[0]
    assert c_dist.wraps_origin is False
    assert ClusterScorer.score_cluster(c_dist).flags == []
    assert ";" in c_dist.format_coordinates(sep=";")


def test_split_mame_nter_cter_opposite_subtypes_not_counted_as_dual_copy():
    """Single split MamE (MamE_Nter + MamE_Cter_greigite) is excluded from both M and dual_copy."""
    # Even with two other regular core genes (MamO, MamP), excluding split MamE keeps M=2 (<3) -> Greigite!
    c = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamA_greigite", "MamA", "core", "Desulfo_greigite"),
        _make_hit("g2", "ctg1", "MamB_greigite", "MamB", "core", "Desulfo_greigite"),
        _make_hit("g3", "ctg1", "MamE_Cter_greigite", "MamE", "core", "Desulfo_greigite"),
        _make_hit("g4", "ctg1", "MamE_Nter", "MamE", "core", "Alpha_classical"),
        _make_hit("g5", "ctg1", "MamO_domain", "MamO", "core", "Alpha_classical"),
        _make_hit("g6", "ctg1", "MamP", "MamP", "core", "Alpha_classical"),
    ])
    sc = ClusterScorer.score_cluster(c)
    assert sc.dual_copy_families == []
    assert sc.magnetite_core_families == ["MamO", "MamP"]
    assert sc.mineral_type == "Greigite"


def test_mineral_type_policy_b_complete_magnetite_with_1_or_2_greigite_hits_stays_magnetite():
    """Policy B: 1 <= G <= 2 with M >= 3 stays Magnetite, while 1 <= G <= 2 with M < 3 is Putative_Greigite."""
    # Complete magnetite MGC (M=4 excl. MamK) + 1 divergent MamI_greigite -> Magnetite
    c_mag = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamA", "MamA", "core", "Alpha_classical"),
        _make_hit("g2", "ctg1", "MamB", "MamB", "core", "Alpha_classical"),
        _make_hit("g3", "ctg1", "MamE_full", "MamE", "core", "Alpha_classical"),
        _make_hit("g4", "ctg1", "MamM", "MamM", "core", "Alpha_classical"),
        _make_hit("g5", "ctg1", "MamI_greigite", "MamI", "core", "Desulfo_greigite"),
    ])
    assert ClusterScorer.score_cluster(c_mag).mineral_type == "Magnetite"

    # Partial greigite MGC (G=2, M=2 excl. MamK, plus shared MamK) -> Putative_Greigite
    c_put = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamB_greigite", "MamB", "core", "Desulfo_greigite"),
        _make_hit("g2", "ctg1", "MamQ_greigite", "MamQ", "core", "Desulfo_greigite"),
        _make_hit("g3", "ctg1", "MamK", "MamK", "core", "Alpha_classical"),
        _make_hit("g4", "ctg1", "MamO_domain", "MamO", "core", "Alpha_classical"),
        _make_hit("g5", "ctg1", "MamP", "MamP", "core", "Alpha_classical"),
    ])
    assert ClusterScorer.score_cluster(c_put).mineral_type == "Putative_Greigite"


def test_circular_3_subclusters_core_span_and_small_plasmid_origin_bridge():
    """Verify 3-subcluster circular coords, origin core_span_bp, exact 25 kb gap, and small plasmid non-merge."""
    engine = ClusterEngine(max_intergenic_genes=5, max_distance_bp=25_000, min_distinct_families=3, flank_genes=3)

    # 1. 4 Mb circular chromosome with origin-bridged pair (1..5 kb & 3995..3999 kb) + independent island at 2 Mb
    starts = [
        1000, 2000, 3000,
        2_000_000, 2_001_000, 2_002_000,
        3_996_000, 3_997_000, 3_998_000,
    ]
    fams = ["MamA", "MamB", "MamE", "MamK", "MamM", "MamO", "MamI", "MamP", "MamQ"]
    genes, hits = [], {}
    for i, (st, fam) in enumerate(zip(starts, fams)):
        genes.append(PredictedGene(f"g{i}", "chr", st, st + 899, 1, i, "M" * 200))
        hits[f"g{i}"] = _make_hit(f"g{i}", "chr", fam, fam, start=st, idx=i)

    c3 = engine.cluster_genome("c3", genes, hits, {"chr": 4_000_000}, {"chr": True})[0]
    assert c3.wraps_origin is True
    assert c3.wraps_origin_compact is False
    assert c3.span_bp < 20_000  # must NOT collapse the 2 Mb internal gap!
    assert c3.core_span_bp <= c3.span_bp
    coords = c3.format_coordinates(sep=";")
    assert "chr:2000000-2002899" in coords
    assert "chr:3996000-3899 (spans origin)" in coords
    assert ClusterScorer.score_cluster(c3).flags == ["Spans_Origin"]

    # 2. Compact 2-subcluster origin-spanning MGC: core_span_bp is accurately wrapped (not ~95 kb)
    s2 = [1000, 2000, 3000, 40_000, 41_000, 96_000, 97_000, 98_000]
    f2 = {0: "MamA", 1: "MamB", 2: "MamE", 5: "MamK", 6: "MamM", 7: "MamO"}
    g2, h2 = [], {}
    for i, st in enumerate(s2):
        g2.append(PredictedGene(f"g{i}", "chr", st, st + 899, 1, i, "M" * 200))
        if i in f2:
            h2[f"g{i}"] = _make_hit(f"g{i}", "chr", f2[i], f2[i], start=st, idx=i)
    c2 = engine.cluster_genome("c2", g2, h2, {"chr": 100_000}, {"chr": True})[0]
    assert c2.wraps_origin is True and c2.wraps_origin_compact is True
    assert c2.core_span_bp == 100_000 - 96_000 + 1 + 3899
    assert c2.core_span_bp <= c2.span_bp

    # 3. Small 30 kb plasmid: 2+2 origin-bridged hits with 6 non-mag genes inside must NOT merge across interior
    sp_starts = [1000 + i * 2500 for i in range(10)]
    sp_fams = {0: "MamA", 1: "MamB", 8: "MamE", 9: "MamQ"}
    sp_g, sp_h = [], {}
    for i, st in enumerate(sp_starts):
        sp_g.append(PredictedGene(f"p{i}", "pl", st, st + 899, 1, i, "M" * 200))
        if i in sp_fams:
            sp_h[f"p{i}"] = _make_hit(f"p{i}", "pl", sp_fams[i], sp_fams[i], start=st, idx=i)
    c_pl = engine.cluster_genome("pl", sp_g, sp_h, {"pl": 26_000}, {"pl": True})[0]
    assert len(c_pl.subclusters) == 2
    assert c_pl.wraps_origin is True
    assert ClusterScorer.score_cluster(c_pl).flags == ["Spans_Origin"]

    # 4. Exact 25,000 bp origin gap rescue consistency
    eng0 = ClusterEngine(max_intergenic_genes=5, max_distance_bp=25_000, min_distinct_families=3, flank_genes=0)
    e_starts = [15_001, 16_001, 88_000, 89_001]  # last_end = 90_000, first_start = 15_001 -> gap = 10000+15000 = 25000
    e_fams = ["MamA", "MamB", "MamE", "MamQ"]
    eg, eh = [], {}
    for i, (st, fam) in enumerate(zip(e_starts, e_fams)):
        en = 90_000 if i == 3 else st + 899
        eg.append(PredictedGene(f"e{i}", "chr", st, en, 1, i, "M" * 200))
        eh[f"e{i}"] = _make_hit(f"e{i}", "chr", fam, fam, start=st, idx=i)
        eh[f"e{i}"].end = en
    c_exact = eng0.cluster_genome("ex", eg, eh, {"chr": 100_000}, {"chr": True})[0]
    assert c_exact.wraps_origin is True


def test_customize_clinker_html_external_note_and_null_strand(tmp_path):
    """External GenBank /note must not block coordinate fallback, and strand=None must not crash."""
    import json
    from magcluster.viz.clinker_export import customize_clinker_html

    hit = _make_hit("ctg1_1", "ctg1", "MamK", "MamK", start=101, idx=0)
    hit.end = 1000
    fake_data = {
        "clusters": [{
            "name": "ext",
            "loci": [{
                "name": "ctg1",
                "genes": [{
                    "uid": "u1",
                    "start": 100,
                    "end": 1000,
                    "strand": None,
                    "names": {"note": "Derived by automated computational analysis"},
                    "label": "",
                }],
            }],
        }],
        "groups": [],
    }
    html = tmp_path / "ext.html"
    html.write_text(f"<script>const data={json.dumps(fake_data)};function serialise(){{}}</script>", encoding="utf-8")
    assert customize_clinker_html(html, {"ctg1_1": hit}) is True
    out = json.loads(re.search(r"<script>const data=(.*?);function serialise", html.read_text()).group(1))
    g = out["clusters"][0]["loci"][0]["genes"][0]
    assert g["is_mag"] is True
    assert g["label"] == "MamK"


def test_logger_detail_preserves_contig_edge_bracket_literal(tmp_path):
    """detail() must keep literal '[contig edge]' and bracketed paths in magcluster.log."""
    from magcluster.utils.logger import attach_log_file, detach_log_file, detail, say

    log_path = tmp_path / "test.log"
    h = attach_log_file(log_path)
    try:
        say("[bold green]Hello[/bold green]")
        detail("  sp_loc1: ctg1:1-900 (3 genes; 3 families) [contig edge]")
        detail("  input: /data/bin[v1].fna")
    finally:
        detach_log_file(h)

    text = log_path.read_text(encoding="utf-8")
    assert "[bold green]" not in text
    assert "Hello" in text
    assert "[contig edge]" in text
    assert "/data/bin[v1].fna" in text


def test_contig_only_genbank_and_external_mgc_title_and_clinker_force_guard(tmp_path):
    """Verify UndefinedSequence GenBank safety, external MGC title HMM search, and clinker --force guard."""
    from Bio import SeqIO
    from magcluster.core.parser import parse_genbank
    from magcluster.core.predict import GenePredictor
    from magcluster.viz.clinker_export import export_cluster_to_genbank, inspect_genbank_for_clinker

    # 1. CONTIG-only GenBank (no ORIGIN) with 0 CDSs -> needs_prediction is False, predictor returns []
    contig_only_gbk = tmp_path / "contig_only.gbff"
    contig_only_gbk.write_text(
        "LOCUS       NZ_TEST                 5000 bp    DNA     linear   BCT 01-JAN-2026\n"
        "DEFINITION  Test contig-only record.\n"
        "ACCESSION   NZ_TEST\n"
        "VERSION     NZ_TEST.1\n"
        "FEATURES             Location/Qualifiers\n"
        "CONTIG      join(AAAA01000001.1:1..5000)\n"
        "//\n",
        encoding="utf-8",
    )
    parsed = parse_genbank(contig_only_gbk)
    assert parsed.needs_prediction is False
    assert GenePredictor().predict_from_records(list(parsed.contigs.values())) == []

    # Exporting a cluster backed by an UndefinedSequence contig safely falls back to N's
    h1 = _make_hit("g1", "NZ_TEST.1", "MamK", "MamK", start=1, idx=0)
    c = _make_cluster_from_hits([h1], genome_name="undef")
    c.subclusters[0].contig_id = "NZ_TEST.1"
    c.subclusters[0].all_genes[0].contig_id = "NZ_TEST.1"
    out_gbk = tmp_path / "undef_MGC.gbk"
    export_cluster_to_genbank(c, parsed.contigs, {"g1": h1}, out_gbk)
    out_rec = SeqIO.read(str(out_gbk), "genbank")
    assert str(out_rec.seq).startswith("NNN")

    # 2. External GenBank whose DEFINITION contains "Magnetosome Gene Cluster" (no MagCluster /note)
    expected_gbk = Path(__file__).resolve().parent / "data" / "expected_1.0.0" / "MSR-1_MAI90kb_MGC.gbk"
    ext_titled = tmp_path / "published_mgc.gbk"
    recs = list(SeqIO.parse(str(expected_gbk), "genbank"))
    for r in recs:
        r.description = "Magnetospirillum sp. Magnetosome Gene Cluster, complete sequence"
        for f in r.features:
            if f.type == "CDS":
                f.qualifiers.pop("note", None)
    SeqIO.write(recs, str(ext_titled), "genbank")
    ext_hits, _, is_mc = inspect_genbank_for_clinker(ext_titled)
    assert is_mc is False
    assert len(ext_hits) >= 25

    # 3. `magcluster clinker --force` with unmatched -l list must NOT delete existing HTML
    plot_dir = tmp_path / "plots"
    plot_dir.mkdir()
    existing = plot_dir / "MGC_align.html"
    existing.write_text("<html>keep me</html>", encoding="utf-8")
    bad_list = tmp_path / "bad.txt"
    bad_list.write_text("NON_EXISTENT_CLUSTER\n", encoding="utf-8")
    res = runner.invoke(app, ["clinker", str(expected_gbk.parent), "-l", str(bad_list), "-o", str(plot_dir), "--force"])
    assert res.exit_code == 1
    assert existing.exists() and existing.read_text() == "<html>keep me</html>"


def test_read_gbk_record_names_prefers_version_and_supports_gzip(tmp_path):
    import gzip
    from magcluster.viz.clinker_export import _read_gbk_record_names

    gbk_text = (
        "LOCUS       NC_023065               1000 bp    DNA     linear   BCT 01-JAN-2026\n"
        "ACCESSION   NC_023065\n"
        "VERSION     NC_023065.1\n"
        "//\n"
    )
    plain = tmp_path / "plain.gbk"
    plain.write_text(gbk_text, encoding="utf-8")
    gz = tmp_path / "comp.gbk.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as f:
        f.write(gbk_text)

    assert _read_gbk_record_names(plain) == [("NC_023065", "NC_023065.1")]
    assert _read_gbk_record_names(gz) == [("NC_023065", "NC_023065.1")]


# --------------------------------------------------------------------------
# 5. Fresh-audit regression tests (circular=false, AMB-1 accession, flank
#    overlap partitioning, pure-greigite Both guards, multi-contig circular,
#    empty translations, corrupt gzip, and API io deduplication)
# --------------------------------------------------------------------------
def test_circular_fasta_header_rejects_circular_false_and_accepts_is_circular_true():
    from magcluster.core.parser import _is_circular_fasta_header

    assert _is_circular_fasta_header("contig_1 length=4000000 circular=true") is True
    assert _is_circular_fasta_header("contig_1 length=4000000 is_circular=true") is True
    assert _is_circular_fasta_header("contig_1 length=4000000 circular=Y") is True
    assert _is_circular_fasta_header("contig_1 length=4000000 circular=1") is True
    assert _is_circular_fasta_header("contig_1 circular chromosome") is True
    assert _is_circular_fasta_header("CP012345.1 Desulfovibrio sp. Rhode Island, complete genome") is True

    assert _is_circular_fasta_header("contig_2 length=50000 depth=1.00x circular=false") is False
    assert _is_circular_fasta_header("contig_2 length=50000 circular=no") is False
    assert _is_circular_fasta_header("contig_2 length=50000 circular=N") is False
    assert _is_circular_fasta_header("contig_2 length=50000 circular=0") is False
    assert _is_circular_fasta_header("contig_2 length=50000 is_circular=false") is False
    assert _is_circular_fasta_header("contig_2 non-circular chromosome") is False
    assert _is_circular_fasta_header("contig_2 complete sequence circular=false") is False
    assert _is_circular_fasta_header("contig_2 circular=unknown") is False


def test_accession_regex_preserves_model_strain_names_like_amb1_and_extracts_real_accessions():
    from magcluster.viz.clinker_export import shorten_gbk_record_ids

    # Ordinary strain/sample names must NOT be misparsed as INSDC accessions
    assert shorten_gbk_record_ids("AMB-1", "AMB-1_MGC", 1, 1) == ("AMB-1_MGC", "AMB-1_MGC")
    assert shorten_gbk_record_ids("AMB-1", "AMB-1_MGC", 2, 4) == ("AMB-1_loc2", "AMB-1_loc2")
    assert shorten_gbk_record_ids("BACTERIUM_01", "BACTERIUM_01_MGC", 1, 1) == ("BACTERIUM_01_MGC", "BACTERIUM_01_MGC")
    assert shorten_gbk_record_ids("FRESHWATER_MAG12", "FRESHWATER_MAG12_MGC", 1, 1) == (
        "FRESHWATER_MAG12_MGC",
        "FRESHWATER_MAG12",
    )
    assert shorten_gbk_record_ids("OUTGROUP_strain3", "OUTGROUP_strain3_MGC", 1, 1)[0] == "OUTGROUP_strain3_MGC"
    assert shorten_gbk_record_ids("OXIC_bin_04", "OXIC_bin_04_MGC", 1, 1)[0] == "OXIC_bin_04_MGC"

    # Real INSDC / RefSeq accessions are still extracted accurately
    assert shorten_gbk_record_ids("GCF_000009985.1_ASM998v1_genomic", "x_MGC", 1, 2)[0] == "GCF_000009985.1_loc1"
    assert shorten_gbk_record_ids("NC_007626.1_AMB-1", "x_MGC", 1, 1)[0] == "NC_007626.1_MGC"
    assert shorten_gbk_record_ids("NZ_CP012345.1_strain", "x_MGC", 1, 1)[0] == "NZ_CP012345.1_MGC"
    assert shorten_gbk_record_ids("CP012345.1_MSR1", "x_MGC", 1, 1)[0] == "CP012345.1_MGC"
    assert shorten_gbk_record_ids("JADFXP010000009.1_bin", "x_MGC", 1, 1)[0] == "JADFXP010000009.1_MGC"


def test_origin_bridged_overlapping_flanks_partitioned_disjointly_without_duplication():
    """On a small circular replicon with large --flank, origin-bridged slices never duplicate genes or exceed contig_len."""
    starts = [1000 + i * 2000 for i in range(10)]
    fams = {0: "MamA", 1: "MamB", 8: "MamE", 9: "MamQ"}
    genes, hits = [], {}
    for i, st in enumerate(starts):
        genes.append(PredictedGene(f"p{i}", "pl", st, st + 899, 1, i, "M" * 200))
        if i in fams:
            hits[f"p{i}"] = _make_hit(f"p{i}", "pl", fams[i], fams[i], start=st, idx=i)

    eng = ClusterEngine(max_intergenic_genes=5, max_distance_bp=25_000, min_distinct_families=3, flank_genes=5)
    c = eng.cluster_genome("pl", genes, hits, {"pl": 22_000}, {"pl": True})[0]
    assert len(c.subclusters) == 2
    ids0 = [g.gene_id for g in c.subclusters[0].all_genes]
    ids1 = [g.gene_id for g in c.subclusters[1].all_genes]
    assert set(ids0).isdisjoint(set(ids1))
    assert len(c.all_genes) == len(set(g.gene_id for g in c.all_genes))
    assert c.wraps_origin is True and c.wraps_origin_compact is True
    assert c.start_bp > c.end_bp
    assert c.span_bp <= c.contig_length
    assert c.core_span_bp <= c.span_bp


def test_multi_contig_circular_chromosome_plus_plasmid_retains_spans_origin():
    """Origin-spanning MGC on circular chr1 plus a second locus on plasmid1 reports both Multi_Contig and Spans_Origin."""
    chr_starts = [500, 1600, 3000, 50_000, 96_000, 97_500, 98_600]
    chr_fams = {0: "MamA", 1: "MamB", 5: "MamE", 6: "MamQ"}
    genes, hits = [], {}
    for i, st in enumerate(chr_starts):
        genes.append(PredictedGene(f"c{i}", "chr1", st, st + 899, 1, i, "M" * 200))
        if i in chr_fams:
            hits[f"c{i}"] = _make_hit(f"c{i}", "chr1", chr_fams[i], chr_fams[i], start=st, idx=i)

    pl_g, pl_h = _build_contig(
        "plasmid1",
        5,
        {
            1: ("Mad1", "Mad1", "accessory", "Desulfo_Mad"),
            2: ("Mad2", "Mad2", "accessory", "Desulfo_Mad"),
            3: ("Mad10", "Mad10", "accessory", "Desulfo_Mad"),
        },
    )
    eng = ClusterEngine(max_intergenic_genes=5, max_distance_bp=25_000, min_distinct_families=3, flank_genes=1)
    c = eng.cluster_genome(
        "sp",
        genes + pl_g,
        {**hits, **pl_h},
        {"chr1": 100_000, "plasmid1": 20_000},
        {"chr1": True, "plasmid1": True},
    )[0]
    assert c.is_multi_contig is True
    assert c.wraps_origin is True
    assert "chr1" in c.origin_contigs
    coords = c.format_coordinates(sep=";")
    assert "chr1:96000-3899 (spans origin)" in coords
    assert "plasmid1:1-4900" in coords
    sc = ClusterScorer.score_cluster(c)
    assert sc.flags == ["Multi_Contig", "Spans_Origin"]
    assert c.span_bp < 20_000


def test_mineral_type_pure_greigite_tandem_mamt_and_single_copy_stay_greigite():
    """Pure greigite MGCs with tandem MamT or single-copy regular hits (no MamL/MamM, no core dual copy) stay Greigite."""
    # Case 1: Tandem MamT + MamT_greigite (accessory) + MamO_domain + MamP (M=2, no core dual copy) -> Greigite
    c_mamt = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamA_greigite", "MamA", "core", "Desulfo_greigite"),
        _make_hit("g2", "ctg1", "MamB_greigite", "MamB", "core", "Desulfo_greigite"),
        _make_hit("g3", "ctg1", "MamQ_greigite", "MamQ", "core", "Desulfo_greigite"),
        _make_hit("g4", "ctg1", "MamT_greigite", "MamT", "accessory", "Desulfo_greigite"),
        _make_hit("g5", "ctg1", "MamT", "MamT", "accessory", "Alpha_classical"),
        _make_hit("g6", "ctg1", "MamO_domain", "MamO", "core", "Alpha_classical"),
        _make_hit("g7", "ctg1", "MamP", "MamP", "core", "Alpha_classical"),
    ])
    sc_mamt = ClusterScorer.score_cluster(c_mamt)
    assert sc_mamt.dual_copy_families == []
    assert sc_mamt.mineral_type == "Greigite"

    # Case 2: Single-copy divergent pure greigite (G=4, M=3 from MamA/MamO/MamP, 0 dual copies, no MamL/MamM) -> Greigite
    c_single = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamB_greigite", "MamB", "core", "Desulfo_greigite"),
        _make_hit("g2", "ctg1", "MamE_greigite", "MamE", "core", "Desulfo_greigite"),
        _make_hit("g3", "ctg1", "MamQ_greigite", "MamQ", "core", "Desulfo_greigite"),
        _make_hit("g4", "ctg1", "MamT_greigite", "MamT", "accessory", "Desulfo_greigite"),
        _make_hit("g5", "ctg1", "MamA", "MamA", "core", "Alpha_classical"),
        _make_hit("g6", "ctg1", "MamO_domain", "MamO", "core", "Alpha_classical"),
        _make_hit("g7", "ctg1", "MamP", "MamP", "core", "Alpha_classical"),
        _make_hit("g8", "ctg1", "MamK", "MamK", "core", "Alpha_classical"),
    ])
    assert ClusterScorer.score_cluster(c_single).mineral_type == "Greigite"

    # Case 3: Disjoint fragmented dual MAG with magnetite-specific MamL/MamM (G=3, M=3 incl. MamL/MamM) -> Both
    c_disjoint_both = _make_cluster_from_hits([
        _make_hit("g1", "ctg1", "MamA_greigite", "MamA", "core", "Desulfo_greigite"),
        _make_hit("g2", "ctg1", "MamB_greigite", "MamB", "core", "Desulfo_greigite"),
        _make_hit("g3", "ctg1", "MamE_greigite", "MamE", "core", "Desulfo_greigite"),
        _make_hit("g4", "ctg2", "MamL", "MamL", "core", "Alpha_classical"),
        _make_hit("g5", "ctg2", "MamM", "MamM", "core", "Alpha_classical"),
        _make_hit("g6", "ctg2", "MamO_full", "MamO", "core", "Alpha_classical"),
    ])
    assert ClusterScorer.score_cluster(c_disjoint_both).mineral_type == "Both"


def test_empty_translations_skipped_and_corrupt_gzip_raises_input_format_error(tmp_path):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    from magcluster.core.parser import InputFormatError, parse_faa, parse_genbank, parse_input

    # 1. Stop-only CDS without /translation produces empty string after rstrip('*') -> skipped
    gbk = tmp_path / "stop_only.gbk"
    rec = SeqRecord(Seq("TAACTG"), id="ctg1", name="ctg1")
    rec.annotations["molecule_type"] = "DNA"
    rec.features = [SeqFeature(FeatureLocation(0, 3, strand=1), type="CDS", qualifiers={"locus_tag": ["STOP1"]})]
    SeqIO.write([rec], str(gbk), "genbank")
    assert parse_genbank(gbk).genes == []

    # 2. Empty record inside a multi-record .faa is skipped cleanly
    faa = tmp_path / "with_empty.faa"
    faa.write_text(">empty_rec\n*\n>valid_rec\nMSEGEGQAKNRL\n", encoding="utf-8")
    parsed_faa = parse_faa(faa)
    assert [g.gene_id for g in parsed_faa.genes] == ["valid_rec"]

    # 3. Corrupted gzip raises InputFormatError (not raw zlib.error)
    bad_gz = tmp_path / "corrupt.fna.gz"
    bad_gz.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03\xff\xff\xff\xff")
    with pytest.raises(InputFormatError):
        parse_input(bad_gz)


def test_io_table_deduplication_and_api_unclustered_export_without_genome_name(tmp_path):
    from magcluster.utils.io import export_genome_magnetosome_fasta, write_genome_cluster_tables

    g1 = PredictedGene("gene_00001", "ctg1", 100, 900, 1, 0, "M" * 200)
    h1 = _make_hit("gene_00001", "ctg1", "MamK", "MamK", start=100, idx=0)

    # 1. Direct Python API call without setting g.genome_name still exports unclustered hits
    tsv = write_genome_cluster_tables("my_strain", [], {"gene_00001": h1}, tmp_path, genes_lookup={"gene_00001": g1})
    faa = export_genome_magnetosome_fasta("my_strain", [], {"gene_00001": h1}, tmp_path, genes_lookup={"gene_00001": g1})
    assert len(tsv.read_text(encoding="utf-8").splitlines()) == 2
    assert ">gene_00001 cluster=unclustered " in faa.read_text(encoding="utf-8")

    # 2. Duplicate gene in cluster.all_genes never writes duplicate TSV rows
    c = _make_cluster_from_hits([h1], genome_name="my_strain")
    c.all_genes = [g1, g1]
    tsv2 = write_genome_cluster_tables("my_strain", [c], {"gene_00001": h1}, tmp_path)
    assert len(tsv2.read_text(encoding="utf-8").splitlines()) == 2


def test_customize_clinker_html_never_leaks_coords_onto_flanking_genes(tmp_path):
    """Non-magnetosome flanking genes with distinct locus_tag must never inherit another contig/locus hit by coords."""
    import json
    from magcluster.viz.clinker_export import customize_clinker_html

    # Trigger A (`magcluster run`): ctg2:1-900 has FeoB (gene_00999); slice loc1 starts with
    # non-magnetosome flanking gene_00050 at slice-local 0..900 (1-based 1..900).
    # Trigger B (`magcluster clinker` multi-slice): loc1 starts with MamA (gene_00001 at 1..900);
    # loc2 starts with non-magnetosome flanking gene_00200 (and even a locus_tag-less hypothetical) at 1..900.
    feo_hit = _make_hit("gene_00999", "ctg2", "FeoB", "FeoB", role="iron_context", lineage="Iron_context", start=1, idx=0)
    feo_hit.end = 900
    mama_hit = _make_hit("gene_00001", "sp_loc1", "MamA", "MamA", start=1, idx=0)
    mama_hit.end = 900

    fake_data = {
        "clusters": [
            {
                "name": "run_cluster",
                "loci": [{
                    "name": "sp_loc1",
                    "genes": [{
                        "uid": "u_flank",
                        "start": 0,
                        "end": 900,
                        "strand": 1,
                        "names": {"locus_tag": "gene_00050", "gene": ""},
                        "label": "",
                    }],
                }],
            },
            {
                "name": "multi_slice",
                "loci": [
                    {
                        "name": "sp_loc1",
                        "genes": [{
                            "uid": "u_mama",
                            "start": 0,
                            "end": 900,
                            "strand": 1,
                            "names": {"locus_tag": "gene_00001", "gene": "MamA"},
                            "label": "MamA",
                        }],
                    },
                    {
                        "name": "sp_loc2",
                        "genes": [
                            {
                                "uid": "u_flank_lt",
                                "start": 0,
                                "end": 900,
                                "strand": 1,
                                "names": {"locus_tag": "gene_00200", "gene": ""},
                                "label": "",
                            },
                            {
                                "uid": "u_flank_nolt",
                                "start": 0,
                                "end": 900,
                                "strand": 1,
                                "names": {"gene": ""},
                                "label": "",
                            },
                        ],
                    },
                ],
            },
        ],
        "groups": [],
    }
    html = tmp_path / "no_leak.html"
    html.write_text(f"<script>const data={json.dumps(fake_data)};function serialise(){{}}</script>", encoding="utf-8")
    assert customize_clinker_html(
        html,
        {},
        hits_by_cluster={
            "run_cluster": {"gene_00999": feo_hit},
            "multi_slice": {"gene_00001": mama_hit},
        },
    ) is True
    out = json.loads(re.search(r"<script>const data=(.*?);function serialise", html.read_text()).group(1))
    g_run = out["clusters"][0]["loci"][0]["genes"][0]
    assert g_run["is_mag"] is False and g_run["label"] == ""

    g_loc1 = out["clusters"][1]["loci"][0]["genes"][0]
    g_loc2_lt, g_loc2_nolt = out["clusters"][1]["loci"][1]["genes"]
    assert g_loc1["is_mag"] is True and g_loc1["label"] == "MamA"
    assert g_loc2_lt["is_mag"] is False and g_loc2_lt["label"] == ""
    assert g_loc2_nolt["is_mag"] is False and g_loc2_nolt["label"] == ""


def test_small_circular_origin_bridged_mgc_preserves_spans_origin_when_interior_shorter():
    """2+2 MGC rescued solely by circular origin bridge keeps wraps_origin=True even if interior_bp < origin_bp."""
    # 20 kb circular plasmid:
    # g0, g1: MamA, MamB at 1,000..2,899
    # g2..g7: 6 compact non-MGC genes (gene_gap=6 > 5, internal_bp_gap=6,100 bp)
    # g8, g9: MamE, MamQ at 9,000..10,899
    # g10: 1 non-MGC gene at 15,000..15,899 (across 0 bp: origin_gene_gap=1 <= 5, origin_bp_gap=10,100 <= 25,000)
    coords = [
        (1000, 1899), (2000, 2899),
        (3000, 3899), (4000, 4899), (5000, 5899), (6000, 6899), (7000, 7899), (8000, 8899),
        (9000, 9899), (10000, 10899),
        (15000, 15899),
    ]
    fams = {0: "MamA", 1: "MamB", 8: "MamE", 9: "MamQ"}
    genes, hits = [], {}
    for i, (st, en) in enumerate(coords):
        genes.append(PredictedGene(f"p{i}", "pl", st, en, 1, i, "M" * 200))
        if i in fams:
            hits[f"p{i}"] = _make_hit(f"p{i}", "pl", fams[i], fams[i], start=st, idx=i)
            hits[f"p{i}"].end = en

    engine = ClusterEngine(max_intergenic_genes=5, max_distance_bp=25_000, min_distinct_families=3, flank_genes=3)
    # Linear: neither side reaches 3 families -> no MGC
    assert engine.cluster_genome("pl", genes, hits, {"pl": 20_000}, {"pl": False}) == []

    # Circular: rescued via origin bridge -> must preserve wraps_origin=True and Spans_Origin!
    res = engine.cluster_genome("pl", genes, hits, {"pl": 20_000}, {"pl": True})
    assert len(res) == 1
    c = res[0]
    assert c.wraps_origin is True
    assert c.wraps_origin_compact is True
    assert c.span_bp <= 20_000
    assert c.core_span_bp <= c.span_bp
    assert "(spans origin)" in c.format_coordinates()
    assert ClusterScorer.score_cluster(c).flags == ["Spans_Origin"]



