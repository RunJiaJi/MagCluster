"""Unit tests for PyHMMER search."""

import pytest
from magcluster.core.predict import PredictedGene
from magcluster.core.search import MagnetosomeSearcher

def test_searcher_initialization():
    searcher = MagnetosomeSearcher()
    assert len(searcher.hmms) > 80
    assert "MamA" in searcher.rules.get("core_families", [])

def test_search_mamk_protein():
    searcher = MagnetosomeSearcher()
    # Canonical MamK sequence from Magnetospirillum gryphiswaldense MSR-1
    mamk_seq = (
        "MSEGEGQAKNRLFLGIDLGTSHTAVMTSRGKKFLLKSVVGYPKDVIGLKLLGRPYVVGDE"
        "AFEMRSYLDLRYPLQDGVLSEISDRDIEVARHLLTHVVKSAEPGANDEICAVIGVPARAS"
        "GANKALLLKMAQEVVHTALVVSEPFMVGYGLDKLNNTIIVDIGAGTTDICALKGTVPGPE"
        "DQVTLTKAGNYLDERLQNAILERHPELQMNTNVACAVKEQFSFVGARGEAATFEFRAAGK"
        "PVRCDVTESVKIACEALMPDIIESIEILLRSFQPEYQATVLQNIVFAGGGSRIRGLAAYV"
        "KDKLRPFGNADVTCVKDPTFDGCRGALRLAEELPPQYWCQLGDVSGQ"
    )
    gene = PredictedGene(
        gene_id="test_mamK",
        contig_id="ctg1",
        start=1,
        end=len(mamk_seq)*3,
        strand=1,
        gene_index=0,
        translation=mamk_seq
    )
    hits = searcher.search_genes([gene])
    assert "test_mamK" in hits
    hit = hits["test_mamK"]
    assert hit.parent_family == "MamK"
    assert hit.role == "core"
    assert hit.score > 200.0
