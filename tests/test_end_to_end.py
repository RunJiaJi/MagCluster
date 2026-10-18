"""End-to-end integration tests on real genomes."""

from pathlib import Path
import pytest
from magcluster.core.cluster import ClusterEngine
from magcluster.core.parser import parse_input
from magcluster.core.scoring import ClusterScorer
from magcluster.core.search import MagnetosomeSearcher

MSR1_GBFF = Path("/Users/sybil/.gemini/antigravity-ide/brain/0bc121af-6698-4630-998e-7a76d2146f06/scratch/controls/MSR-1_genome.gbff")

@pytest.mark.skipif(not MSR1_GBFF.exists(), reason="MSR-1 control genome file not found")
def test_msr1_end_to_end():
    # 1. Parse GBFF
    parsed = parse_input(MSR1_GBFF)
    assert len(parsed.genes) > 4000
    assert not parsed.needs_prediction

    # 2. Search
    searcher = MagnetosomeSearcher()
    hits = searcher.search_genes(parsed.genes)
    assert len(hits) >= 40

    # 3. Cluster
    clusterer = ClusterEngine()
    clusters = clusterer.cluster_genome("MSR-1", parsed.genes, hits)
    assert len(clusters) == 1

    c = clusters[0]
    assert c.span_bp > 20_000
    assert len(c.distinct_core_families) == 10

    # 4. Score
    score = ClusterScorer.score_cluster(c)
    assert score.completeness_score == 1.0
    assert score.confidence == "High"
    assert score.mineral_type == "Magnetite"
    assert score.lineage_signature == "Alphaproteobacteria"
