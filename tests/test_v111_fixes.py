"""Tests for MagCluster v1.1.1 bug fixes and UX improvements."""

import json
import re
import shutil
import sys
import warnings
from pathlib import Path

import pytest
from Bio import BiopythonWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from typer.testing import CliRunner

from magcluster.cli import app
from magcluster.core.cluster import ClusterEngine, MGCCluster, MGCSubcluster
from magcluster.core.parser import parse_input
from magcluster.core.predict import PredictedGene
from magcluster.core.search import MagnetosomeHit
from magcluster.viz.clinker_export import export_cluster_to_genbank, shorten_gbk_record_ids

DATA = Path(__file__).resolve().parent / "data"
MAI_GZ = DATA / "MSR-1_MAI90kb.fna.gz"
EXPECTED_GBK = DATA / "expected_1.0.0" / "MSR-1_MAI90kb_MGC.gbk"

runner = CliRunner()
HAS_CLINKER = shutil.which("clinker") is not None or (Path(sys.executable).parent / "clinker").exists()

ANSI_RE = re.compile(r"\x1b\[[0-9;?]*[ -/]*[@-~]")

MAMK_SEQ = (
    "MSEGEGQAKNRLFLGIDLGTSHTAVMTSRGKKFLLKSVVGYPKDVIGLKLLGRPYVVGDE"
    "AFEMRSYLDLRYPLQDGVLSEISDRDIEVARHLLTHVVKSAEPGANDEICAVIGVPARAS"
    "GANKALLLKMAQEVVHTALVVSEPFMVGYGLDKLNNTIIVDIGAGTTDICALKGTVPGPE"
    "DQVTLTKAGNYLDERLQNAILERHPELQMNTNVACAVKEQFSFVGARGEAATFEFRAAGK"
    "PVRCDVTESVKIACEALMPDIIESIEILLRSFQPEYQATVLQNIVFAGGGSRIRGLAAYV"
    "KDKLRPFGNADVTCVKDPTFDGCRGALRLAEELPPQYWCQLGDVSGQ"
)

FEOB_SEQ = (
    "MSTIKVAVSGQPNCGKSTMFNAITGGSARVGNYPGITVDRLEGVYQQNGDKIHLVDLPGT"
    "YSLTSYSMEELVARNVIVDEHPDVVINMIDATALERSLYLAVQFMEIGAPMVLGLNMMDE"
    "VKRNGVTIDVKKLSTLLGVPVVPCVARVGQGREELMKAVAAAYKASNGKWKPIELSYGPE"
    "LDPVIERMTAMIAEAGFLTGRHDPRWVAIKYLENDEQVIEEGRAAGSLSAKLEALCKEAE"
)


def flat(text: str) -> str:
    return " ".join(ANSI_RE.sub("", text).split())


def run_cli(*args):
    res = runner.invoke(app, [str(a) for a in args])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert "Traceback" not in res.output
    return res


# --------------------------------------------------------------------------
# 1. Protein FASTA named .fasta / .fa auto-detected & duplicate IDs resolved
# --------------------------------------------------------------------------
def test_protein_fasta_with_fasta_extension_and_duplicate_pipe_headers(tmp_path):
    prot_fasta = tmp_path / "representative_magProteins.fasta"
    prot_fasta.write_text(
        f">genbankID| Magnetosome protein MamK AGG16218.1 [FH-1]\n{MAMK_SEQ}\n"
        f">genbankID| Magnetosome protein FeoB1 AGG16234.1 [FH-1]\n{FEOB_SEQ}\n"
        f">genbankID| Magnetosome protein MamK2 AGG16219.1 [FH-1]\n{MAMK_SEQ}\n",
        encoding="utf-8",
    )

    parsed = parse_input(prot_fasta)
    assert parsed.file_type == "faa"
    assert not parsed.needs_prediction
    assert len(parsed.genes) == 3
    assert [g.gene_id for g in parsed.genes] == [
        "genbankID|AGG16218.1",
        "genbankID|AGG16234.1",
        "genbankID|AGG16219.1",
    ]
    assert any("amino acid sequences" in n for n in parsed.notices)
    assert any("duplicate protein sequence ID" in n for n in parsed.notices)

    out_dir = tmp_path / "out"
    res = run_cli("run", prot_fasta, "-o", out_dir, "--no-plot")
    assert res.exit_code == 0
    txt = flat(res.output)
    assert "amino acid sequences" in txt
    tsv = out_dir / "representative_magProteins_mgc_genes.tsv"
    assert tsv.exists()
    rows = tsv.read_text(encoding="utf-8").splitlines()[1:]
    # All 3 magnetosome proteins are retained rather than collapsing into a single 'genbankID|' key
    assert len(rows) == 3
    assert {r.split("\t")[2] for r in rows} == {
        "genbankID|AGG16218.1",
        "genbankID|AGG16234.1",
        "genbankID|AGG16219.1",
    }


def test_duplicate_contig_ids_in_nucleotide_fasta_resolved(tmp_path):
    fna = tmp_path / "dup_ctg.fna"
    dna = "ATG" + "GCT" * 80 + "TAA"
    fna.write_text(f">contigA first\n{dna}\n>contigA second\n{dna}\n", encoding="utf-8")
    parsed = parse_input(fna)
    assert parsed.file_type == "fasta"
    assert list(parsed.contigs.keys()) == ["contigA", "contigA_2"]
    assert any("duplicate contig ID" in n for n in parsed.notices)


# --------------------------------------------------------------------------
# 2. GenBank partial codons & pseudogenes emit zero BiopythonWarning
# --------------------------------------------------------------------------
def test_genbank_partial_codon_and_pseudo_cds_emit_no_warning(tmp_path):
    seq = Seq("ATGCGATCGATCGATCGATCGATCGAA" + "AC")  # 29 bp (not divisible by 3)
    rec = SeqRecord(seq, id="NC_TEST.1", name="NC_TEST", description="Test contig")
    rec.annotations["molecule_type"] = "DNA"

    # Pseudogene without /translation and non-multiple-of-3 length
    pseudo_cds = SeqFeature(
        FeatureLocation(0, 29, strand=1),
        type="CDS",
        qualifiers={"locus_tag": ["LT_PSEUDO"], "pseudo": [""]},
    )
    # Truncated edge CDS without /translation and non-multiple-of-3 length + duplicate locus_tag check
    partial_cds_1 = SeqFeature(
        FeatureLocation(0, 29, strand=1),
        type="CDS",
        qualifiers={"locus_tag": ["LT_DUP"], "codon_start": ["2"]},
    )
    partial_cds_2 = SeqFeature(
        FeatureLocation(0, 28, strand=1),
        type="CDS",
        qualifiers={"locus_tag": ["LT_DUP"]},
    )
    rec.features = [pseudo_cds, partial_cds_1, partial_cds_2]

    gbk_path = tmp_path / "partial.gbff"
    SeqIO.write([rec], str(gbk_path), "genbank")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", BiopythonWarning)
        parsed = parse_input(gbk_path)

    assert len(caught) == 0, [str(w.message) for w in caught]
    # Pseudogene without translation was skipped; two duplicate LT_DUP CDSs were both kept cleanly
    assert [g.gene_id for g in parsed.genes] == ["LT_DUP", "LT_DUP_2"]


# --------------------------------------------------------------------------
# 3. Same-contig multi-subcluster coordinates & genome-wide _loc numbering
# --------------------------------------------------------------------------
def test_same_contig_multi_subcluster_coordinates_and_multi_contig_loc_ids():
    # Part A: Same chromosome with 2 distant subclusters (like AMB-1)
    g1 = PredictedGene("g1", "NC_007626.1", 421240, 436056, 1, 0, "M" * 100)
    g2 = PredictedGene("g2", "NC_007626.1", 1010931, 1045256, 1, 500, "M" * 100)
    s1 = MGCSubcluster("amb_loc1", "NC_007626.1", 421240, 436056, 421240, 436056, all_genes=[g1])
    s2 = MGCSubcluster("amb_loc2", "NC_007626.1", 1010931, 1045256, 1010931, 1045256, all_genes=[g2])
    c_same = MGCCluster(
        cluster_id="amb_MGC",
        genome_name="amb",
        contig_id="NC_007626.1",
        start_bp=421240,
        end_bp=1045256,
        cluster_core_start=421240,
        cluster_core_end=1045256,
        all_genes=[g1, g2],
        subclusters=[s1, s2],
    )
    assert c_same.format_coordinates(sep="; ") == "NC_007626.1:421240-436056; NC_007626.1:1010931-1045256"

    # Part B: Multi-contig genome gets sequential _loc1, _loc2 across contigs (not _loc1, _loc1)
    genes = []
    hits = {}
    for c_idx, ctg in enumerate(("ctgA", "ctgB"), 1):
        for f_idx, fam in enumerate(("MamA", "MamB", "MamK")):
            gid = f"{ctg}_{f_idx}"
            st = 100 + f_idx * 1000
            genes.append(PredictedGene(gid, ctg, st, st + 600, 1, f_idx, "M" * 200))
            hits[gid] = MagnetosomeHit(
                gid, ctg, fam, fam, 200.0, 1e-40, "core", "Alpha", "#EA4335", st, st + 600, 1, f_idx, "M" * 200
            )
    clusters = ClusterEngine(flank_genes=0).cluster_genome("magX", genes, hits)
    assert len(clusters) == 1
    assert [sub.subcluster_id for sub in clusters[0].subclusters] == ["magX_loc1", "magX_loc2"]


# --------------------------------------------------------------------------
# 4. Accession-based GBK headers on medium names, clean DEFINITION, untruncated HTML loci
# --------------------------------------------------------------------------
def test_short_accession_genome_uses_accession_and_clean_single_slice_definition(tmp_path):
    gname = "GCF_000009985.1_ASM998v1_genomic"
    rec_id_1, rec_name_1 = shorten_gbk_record_ids(gname, f"{gname}_MGC", slice_idx=1, n_slices=4)
    rec_id_2, rec_name_2 = shorten_gbk_record_ids(gname, f"{gname}_MGC", slice_idx=2, n_slices=4)
    assert (rec_id_1, rec_id_2) == ("GCF_000009985.1_loc1", "GCF_000009985.1_loc2")
    assert rec_name_1 != rec_name_2
    assert rec_name_1.endswith("_loc1") and rec_name_2.endswith("_loc2")

    # Single-slice DEFINITION has no duplicate "Cluster cluster"
    g1 = PredictedGene("g1", "JADFXP010000009.1", 100, 999, 1, 0, "M" * 300)
    c1 = MGCCluster(
        cluster_id="GCA_015233465.1_ASM1523346v1_genomic_MGC",
        genome_name="GCA_015233465.1_ASM1523346v1_genomic",
        contig_id="JADFXP010000009.1",
        start_bp=100,
        end_bp=999,
        cluster_core_start=100,
        cluster_core_end=999,
        all_genes=[g1],
    )
    gbk_out = tmp_path / "single.gbk"
    export_cluster_to_genbank(c1, {}, {}, gbk_out)
    rec = SeqIO.read(str(gbk_out), "genbank")
    assert rec.id == "GCA_015233465.1_MGC"
    assert "Cluster cluster" not in rec.description
    assert rec.description.startswith("Magnetosome Gene Cluster from GCA_015233465.1_ASM1523346v1_genomic")


# --------------------------------------------------------------------------
# 5. Batch progress output, MGC_align_batches.tsv manifest, and clinker UX warnings
# --------------------------------------------------------------------------
@pytest.mark.skipif(not HAS_CLINKER, reason="clinker not installed")
def test_batch_progress_not_swallowed_and_batch_index_written_and_cleaned(tmp_path):
    gbk_dir = tmp_path / "gbks"
    gbk_dir.mkdir()
    for tag in ("strain01_MGC.gbk", "strain02_MGC.gbk", "strain03_MGC.gbk"):
        shutil.copy(EXPECTED_GBK, gbk_dir / tag)

    out_dir = tmp_path / "plots"
    res = run_cli("clinker", gbk_dir, "-o", out_dir, "-b", "2", "--clinker-opts", "-na -ufo")
    assert res.exit_code == 0
    out_txt = flat(res.output)
    # Batch label is visible (not stripped by Rich markup)
    assert "(batch 1/2, 2 clusters: strain01_MGC .. strain02_MGC)" in out_txt
    assert "(batch 2/2, 1 clusters: strain03_MGC)" in out_txt

    idx_tsv = out_dir / "MGC_align_batches.tsv"
    assert idx_tsv.exists()
    lines = idx_tsv.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "batch\thtml_file\torder_in_batch\tcluster_id\tgbk_file"
    assert len(lines) == 4
    assert lines[1].startswith("1\tMGC_align_batch01.html\t1\tstrain01_MGC\t")
    assert lines[3].startswith("2\tMGC_align_batch02.html\t1\tstrain03_MGC\t")

    # Verify untruncated locus names in Clinker HTML JSON (MSR-1_MAI90kb_loc1..4 instead of MSR-1_MAI90kb_lo)
    html1 = (out_dir / "MGC_align_batch01.html").read_text(encoding="utf-8")
    d1 = json.loads(re.search(r"<script>const data=(.*?);function serialise", html1).group(1))
    loci_names = [loc["name"] for loc in d1["clusters"][0]["loci"]]
    assert loci_names == [
        "MSR-1_MAI90kb_loc1",
        "MSR-1_MAI90kb_loc2",
        "MSR-1_MAI90kb_loc3",
        "MSR-1_MAI90kb_loc4",
    ]

    # Re-running with -b 0 --force cleans up stale batch HTMLs AND stale MGC_align_batches.tsv
    res_single = run_cli("clinker", gbk_dir, "-o", out_dir, "-b", "0", "--clinker-opts", "-na -ufo", "--force")
    assert res_single.exit_code == 0
    assert (out_dir / "MGC_align.html").exists()
    assert not idx_tsv.exists()
    assert not (out_dir / "MGC_align_batch01.html").exists()


def test_clinker_bare_list_without_directory_gives_actionable_hint(tmp_path):
    lst = tmp_path / "order.txt"
    lst.write_text("MSR-1_genome_MGC\n", encoding="utf-8")
    res = run_cli("clinker", "-l", lst, "-o", tmp_path / "out")
    assert res.exit_code == 1
    assert "pass a GenBank/results directory along with -l" in flat(res.output)


@pytest.mark.skipif(not HAS_CLINKER, reason="clinker not installed")
def test_clinker_warns_when_given_unsliced_whole_genome_genbank(tmp_path):
    # Create a synthetic external GenBank with 505 tiny CDSs and no MagCluster /note
    dna = Seq("ATG" + "GCT" * 5 + "TAA") * 505
    rec = SeqRecord(dna, id="NC_BIG.1", name="NC_BIG", description="Whole bacterial chromosome")
    rec.annotations["molecule_type"] = "DNA"
    rec.features = [
        SeqFeature(
            FeatureLocation(i * 21, (i + 1) * 21, strand=1),
            type="CDS",
            qualifiers={"locus_tag": [f"BIG_{i:04d}"], "translation": ["MAAAAA"]},
        )
        for i in range(505)
    ]
    big_gbk = tmp_path / "whole_genome.gbff"
    SeqIO.write([rec], str(big_gbk), "genbank")

    out_html = tmp_path / "big.html"
    res = run_cli("clinker", big_gbk, "-o", out_html, "--clinker-opts", "-na")
    assert res.exit_code == 0
    assert "appears to be an unsliced whole genome" in flat(res.output)
    assert out_html.exists()
    # Direct .html output target does not drop a stray magcluster_clinker.log without --log
    assert not (tmp_path / "magcluster_clinker.log").exists()
    assert not Path("magcluster_clinker.log").exists()
