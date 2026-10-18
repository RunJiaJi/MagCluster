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

def test_clinker_profile_labels(tmp_path):
    from magcluster.core.cluster import MGCCluster
    from magcluster.core.predict import PredictedGene
    from magcluster.core.search import MagnetosomeHit
    from magcluster.viz.clinker_export import export_cluster_to_genbank, customize_clinker_html
    import json
    from Bio import SeqIO

    gene1 = PredictedGene(gene_id="g_mame", contig_id="ctg1", start=100, end=700, strand=1, gene_index=0, translation="MEFULL")
    cluster = MGCCluster(
        cluster_id="test_cluster",
        genome_name="test_genome",
        contig_id="ctg1",
        start_bp=100,
        end_bp=700,
        cluster_core_start=100,
        cluster_core_end=700,
        all_genes=[gene1]
    )
    hits = {
        "g_mame": MagnetosomeHit(
            gene_id="g_mame", contig_id="ctg1", profile="MamE_Cter", parent_family="MamE",
            score=500.0, evalue=1e-100, role="core", lineage="Alpha", color="#EA4335",
            start=100, end=700, strand=1, gene_index=0, translation="MEFULL"
        )
    }

    # 1. GenBank export
    gbk_file = tmp_path / "test.gbk"
    export_cluster_to_genbank(cluster, {}, hits, gbk_file)
    rec = SeqIO.read(str(gbk_file), "genbank")
    cds_feat = [f for f in rec.features if f.type == "CDS"][0]
    assert cds_feat.qualifiers["label"] == ["MamE_Cter"]
    assert cds_feat.qualifiers["gene"] == ["MamE_Cter"]
    assert cds_feat.qualifiers["product"] == ["MamE_Cter"]
    assert "Profile=MamE_Cter" in cds_feat.qualifiers["note"][0]
    assert "Family=MamE" in cds_feat.qualifiers["note"][0]

    # 2. HTML customization
    fake_data = {
        "clusters": [
            {
                "name": "test_cluster",
                "loci": [
                    {
                        "name": "ctg1",
                        "genes": [
                            {
                                "uid": "gene_uid_1",
                                "names": {"locus_tag": "g_mame", "gene": "MamE_Cter"},
                                "label": "MamE_Cter"
                            }
                        ]
                    }
                ]
            }
        ],
        "groups": []
    }
    html_file = tmp_path / "test.html"
    html_file.write_text(f"<script>const data={json.dumps(fake_data)};function serialise(){{}}</script>", encoding="utf-8")
    
    assert customize_clinker_html(html_file, hits) is True
    updated_html = html_file.read_text(encoding="utf-8")
    import re
    m = re.search(r"<script>const data=(.*?);function serialise", updated_html)
    data = json.loads(m.group(1))
    g = data["clusters"][0]["loci"][0]["genes"][0]
    assert g["label"] == "MamE_Cter"
    assert g["names"]["gene"] == "MamE_Cter"
    assert g["names"]["product"] == "MamE_Cter"
    # Verify group is mam family (red)
    mam_group = [grp for grp in data["groups"] if grp["label"] == "mam"][0]
    assert "gene_uid_1" in mam_group["genes"]

    # 3. Verify _full suffix omission for full-length proteins
    from magcluster.viz.clinker_export import get_display_label
    hit_full = MagnetosomeHit(
        gene_id="g_full", contig_id="ctg1", profile="MamE_full", parent_family="MamE",
        score=900.0, evalue=0.0, role="core", lineage="Alpha", color="#EA4335",
        start=100, end=700, strand=1, gene_index=0, translation="MEFULL"
    )
    hit_o_full = MagnetosomeHit(
        gene_id="g_ofull", contig_id="ctg1", profile="MamO_full", parent_family="MamO",
        score=900.0, evalue=0.0, role="core", lineage="Alpha", color="#EA4335",
        start=100, end=700, strand=1, gene_index=0, translation="MOFULL"
    )
    assert get_display_label(hit_full) == "MamE"
    assert get_display_label(hit_o_full) == "MamO"
    assert get_display_label(hits["g_mame"]) == "MamE_Cter"

def test_force_overwrite_existing_outdir(tmp_path):
    from typer.testing import CliRunner
    from magcluster.cli import app

    runner = CliRunner()
    test_faa = tmp_path / "test.faa"
    test_faa.write_text(">p1\nMSQIKIVGVGPGGIAALIGEARALGISDVTVVDLEEAAGRAVAEA\n", encoding="utf-8")

    out_dir = tmp_path / "results_dir"
    out_dir.mkdir(parents=True)
    # Place an existing file in out_dir
    (out_dir / "existing_file.txt").write_text("dummy", encoding="utf-8")

    # 1. Run without --force -> should fail and instruct user to use --force
    res = runner.invoke(app, ["run", str(test_faa), "-o", str(out_dir), "--no-plot"])
    assert res.exit_code == 1
    assert "already exists" in res.output
    assert "--force" in res.output

    # 2. Run with --force -> should proceed
    res_force = runner.invoke(app, ["run", str(test_faa), "-o", str(out_dir), "--no-plot", "--force"])
    assert res_force.exit_code == 0


