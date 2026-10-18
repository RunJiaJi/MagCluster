"""Unit tests for Pyrodigal gene prediction."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from magcluster.core.predict import GenePredictor

def test_gene_predictor_basic():
    predictor = GenePredictor(meta=True)
    # Complete valid coding sequence with shine-dalgarno and canonical stop
    coding_dna = (
        "AGGAGG" + "TTTT" +
        "ATGAAGAAAATCGCTCTTTTAAGCGCCCTCGCCCTTTCCGCTGCTCTCGCCACCGCTACC"
        "GCCAACGCTGACGAAGTACAGATCGTTCTTGACGGCAAGACCCTCGCCAACGGCAAGGTC"
        "ATCGCCGACGGCACCATCGACATCGCCGCCAACGGCAAGGTTCTGGTCGACAAGTAA"
    )

    record = SeqRecord(Seq(coding_dna), id="test_contig_1")
    genes = predictor.predict_from_records([record], genome_name="test_genome")

    assert len(genes) >= 1
    g = genes[0]
    assert g.contig_id == "test_contig_1"
    assert g.length_aa >= 30
    assert g.gene_id.startswith("test_genome_")
