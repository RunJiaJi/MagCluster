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
    # Each genome yields exactly one consolidated MGC
    assert len(clusters) == 1

    c = clusters[0]
    assert c.cluster_id == "MSR-1_MGC"
    assert c.span_bp > 50_000
    assert len(c.distinct_magnetosome_families) >= 20

    # 4. Score
    score = ClusterScorer.score_cluster(c)
    assert score.mineral_type == "Magnetite"
    assert score.n_mag_families >= 20

def test_protein_fasta_without_coords(tmp_path):
    faa_content = (
        ">protein_1 some header without coords\n"
        "MSQIKIVGVGPGGIAALIGEARALGISDVTVVDLEEAAGRAVAEA\n"
        ">protein_2 another header\n"
        "MKSPVTVALVGNPNAGSTSLFNALTGAQSAVANYQRVTTSLNVRE\n"
    )
    test_faa = tmp_path / "test.faa"
    test_faa.write_text(faa_content, encoding="utf-8")

    parsed = parse_input(test_faa)
    assert parsed.file_type == "faa"
    assert not parsed.has_coordinates
    assert len(parsed.genes) == 2
    assert parsed.genes[0].contig_id == "unplaced"
    assert parsed.genes[0].start == 0

def test_per_genome_exports(tmp_path):
    from magcluster.utils.io import write_all_genome_tables, export_all_genome_fastas
    from magcluster.core.search import MagnetosomeHit
    from magcluster.core.predict import PredictedGene

    gene1 = PredictedGene(gene_id="g1", contig_id="c1", start=100, end=500, strand=1, gene_index=0, translation="MTEST")
    gene1.genome_name = "genomeA"
    gene2 = PredictedGene(gene_id="g2", contig_id="c2", start=200, end=600, strand=-1, gene_index=0, translation="MTEST2")
    gene2.genome_name = "genomeB"

    hits = {
        "g1": MagnetosomeHit(
            gene_id="g1", contig_id="c1", profile="MamA", parent_family="MamA",
            score=100.0, evalue=1e-20, role="core", lineage="Alpha", color="#E63946",
            start=100, end=500, strand=1, gene_index=0, translation="MTEST"
        ),
        "g2": MagnetosomeHit(
            gene_id="g2", contig_id="c2", profile="MamB", parent_family="MamB",
            score=80.0, evalue=1e-15, role="core", lineage="Alpha", color="#E63946",
            start=200, end=600, strand=-1, gene_index=0, translation="MTEST2"
        ),
    }
    lookup = {"g1": gene1, "g2": gene2}

    tsvs = write_all_genome_tables(["genomeA", "genomeB"], [], {}, hits, tmp_path, genes_lookup=lookup)
    faas = export_all_genome_fastas(["genomeA", "genomeB"], [], hits, tmp_path, genes_lookup=lookup)

    assert len(tsvs) == 2
    assert (tmp_path / "genomeA_mgc_genes.tsv").exists()
    assert (tmp_path / "genomeB_mgc_genes.tsv").exists()
    assert len(faas) == 2
    assert (tmp_path / "genomeA_magnetosome_proteins.faa").exists()
    assert (tmp_path / "genomeB_magnetosome_proteins.faa").exists()
