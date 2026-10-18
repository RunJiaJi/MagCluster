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
    # With min_distinct_families=3, MSR-1 yields 4 operons: mamXY, mamMRV, mamAB, and mamGFDC-mms6
    assert len(clusters) == 4

    # The main mamAB cluster
    c_main = [c for c in clusters if len(c.distinct_magnetosome_families) >= 10][0]
    assert c_main.span_bp > 20_000
    assert len(c_main.distinct_magnetosome_families) >= 10

    # 4. Score
    score = ClusterScorer.score_cluster(c_main)
    assert score.mineral_type == "Magnetite"
    assert score.n_mag_families >= 10
