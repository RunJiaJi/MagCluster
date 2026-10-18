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
