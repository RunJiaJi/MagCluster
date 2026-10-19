"""Unit tests for ClusterEngine and ClusterScorer."""

import pytest
from magcluster.core.cluster import ClusterEngine
from magcluster.core.predict import PredictedGene
from magcluster.core.scoring import ClusterScorer
from magcluster.core.search import MagnetosomeHit

def test_cluster_engine_synteny():
    # Build 10 synthetic genes:
    # Gene 0: non-MGC
    # Gene 1: MamA (core)
    # Gene 2: MamB (core)
    # Gene 3: non-MGC
    # Gene 4: MamE (core)
    # Gene 5: MamK (core)
    # Gene 6..9: non-MGC
    genes = []
    for i in range(10):
        genes.append(PredictedGene(
            gene_id=f"g_{i}",
            contig_id="ctg_1",
            start=1 + i * 1000,
            end=900 + i * 1000,
            strand=1,
            gene_index=i,
            translation="M" * 100
        ))

    hits = {
        "g_1": MagnetosomeHit("g_1", "ctg_1", "MamA", "MamA", 100.0, 1e-20, "core", "Alpha_classical", "#E63946", 1001, 1900, 1, 1, "M"*100),
        "g_2": MagnetosomeHit("g_2", "ctg_1", "MamB", "MamB", 100.0, 1e-20, "core", "Alpha_classical", "#E63946", 2001, 2900, 1, 2, "M"*100),
        "g_4": MagnetosomeHit("g_4", "ctg_1", "MamE", "MamE", 100.0, 1e-20, "core", "Alpha_classical", "#E63946", 4001, 4900, 1, 4, "M"*100),
        "g_5": MagnetosomeHit("g_5", "ctg_1", "MamK", "MamK", 100.0, 1e-20, "core", "Alpha_classical", "#E63946", 5001, 5900, 1, 5, "M"*100),
    }

    engine = ClusterEngine(max_intergenic_genes=2, min_distinct_families=3, flank_genes=1)
    clusters = engine.cluster_genome("test_sp", genes, hits)

    assert len(clusters) == 1
    c = clusters[0]
    assert c.cluster_id == "test_sp_MGC"
    assert len(c.distinct_magnetosome_families) == 4
    assert "MamA" in c.distinct_magnetosome_families
    assert "MamK" in c.distinct_magnetosome_families

    score = ClusterScorer.score_cluster(c)
    assert score.mineral_type == "Magnetite"
    assert score.n_mag_families == 4

def test_cluster_engine_subcluster_slicing():
    # Construct 40 genes:
    # Subcluster 1: genes 5, 6, 7 (MamA, MamB, MamE)
    # 20 non-MGC genes: 8..27
    # Subcluster 2: genes 28, 29, 30 (MamK, MamM, MamO)
    genes = []
    for i in range(40):
        genes.append(PredictedGene(
            gene_id=f"g_{i}",
            contig_id="chr1",
            start=1 + i * 1000,
            end=900 + i * 1000,
            strand=1,
            gene_index=i,
            translation="M" * 100
        ))

    hits = {
        "g_5": MagnetosomeHit("g_5", "chr1", "MamA", "MamA", 100.0, 1e-20, "core", "Alpha_classical", "#E63946", 5001, 5900, 1, 5, "M"*100),
        "g_6": MagnetosomeHit("g_6", "chr1", "MamB", "MamB", 100.0, 1e-20, "core", "Alpha_classical", "#E63946", 6001, 6900, 1, 6, "M"*100),
        "g_7": MagnetosomeHit("g_7", "chr1", "MamE", "MamE", 100.0, 1e-20, "core", "Alpha_classical", "#E63946", 7001, 7900, 1, 7, "M"*100),
        "g_28": MagnetosomeHit("g_28", "chr1", "MamK", "MamK", 100.0, 1e-20, "core", "Alpha_classical", "#E63946", 28001, 28900, 1, 28, "M"*100),
        "g_29": MagnetosomeHit("g_29", "chr1", "MamM", "MamM", 100.0, 1e-20, "core", "Alpha_classical", "#E63946", 29001, 29900, 1, 29, "M"*100),
        "g_30": MagnetosomeHit("g_30", "chr1", "MamO", "MamO", 100.0, 1e-20, "core", "Alpha_classical", "#E63946", 30001, 30900, 1, 30, "M"*100),
    }

    engine = ClusterEngine(max_intergenic_genes=5, min_distinct_families=3, flank_genes=2)
    clusters = engine.cluster_genome("test_amb1", genes, hits)

    assert len(clusters) == 1
    c = clusters[0]
    assert len(c.subclusters) == 2

    # Subcluster 1 has genes 3..9 (3 hits + 2 flank on each side = 7 genes)
    sub1 = c.subclusters[0]
    assert [g.gene_id for g in sub1.all_genes] == [f"g_{i}" for i in range(3, 10)]

    # Subcluster 2 has genes 26..32 (3 hits + 2 flank on each side = 7 genes)
    sub2 = c.subclusters[1]
    assert [g.gene_id for g in sub2.all_genes] == [f"g_{i}" for i in range(26, 33)]

    # The 16 genes in between (g_10 to g_25) must NOT be in c.all_genes!
    assert len(c.all_genes) == 14
    for mid_i in range(10, 26):
        assert f"g_{mid_i}" not in [g.gene_id for g in c.all_genes]

