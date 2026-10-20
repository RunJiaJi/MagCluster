"""Regression tests for the bugs fixed in MagCluster 1.0.1.

Most tests run the real CLI on a 90 kb slice of the Magnetospirillum
gryphiswaldense MSR-1 genome that contains the magnetosome island
(``tests/data/MSR-1_MAI90kb.fna.gz``; from NCBI RefSeq NC_023065.1).
"""

import gzip
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from typer.testing import CliRunner

from magcluster.cli import app
from magcluster.core.cluster import ClusterEngine
from magcluster.core.parser import (
    InputFormatError,
    genome_name_from_path,
    is_genome_file,
    parse_input,
    sanitize_genome_name,
)
from magcluster.core.predict import PredictedGene
from magcluster.core.search import MagnetosomeHit, MagnetosomeSearcher

DATA = Path(__file__).resolve().parent / "data"
MAI_GZ = DATA / "MSR-1_MAI90kb.fna.gz"
EXPECTED = DATA / "expected_1.0.0"
REPO = Path(__file__).resolve().parent.parent

runner = CliRunner()
HAS_CLINKER = shutil.which("clinker") is not None or (Path(sys.executable).parent / "clinker").exists()


ANSI_RE = re.compile(r"\x1b\[[0-9;?]*[ -/]*[@-~]")


def strip_ansi(text: str) -> str:
    """Remove ANSI colour/style codes (Rich adds them when FORCE_COLOR or GITHUB_ACTIONS is set)."""
    return ANSI_RE.sub("", text)


def flat(text: str) -> str:
    """Strip ANSI codes and collapse Rich line wrapping, so substring checks do not
    depend on terminal width or colour settings."""
    return " ".join(strip_ansi(text).split())


def run_cli(*args):
    res = runner.invoke(app, [str(a) for a in args])
    # Only clean exits are acceptable: never an uncaught Python exception
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert "Traceback" not in res.output
    return res


@pytest.fixture
def mai_fna(tmp_path):
    """Plain FASTA copy of the MSR-1 MAI slice."""
    p = tmp_path / "in" / "MSR-1_MAI90kb.fna"
    p.parent.mkdir()
    with gzip.open(MAI_GZ, "rt") as src:
        p.write_text(src.read())
    return p


def tsv_rows(path: Path):
    return path.read_text().splitlines()[1:]


# --------------------------------------------------------------------------
# Detection outputs are unchanged versus 1.0.0
# --------------------------------------------------------------------------
def test_outputs_byte_identical_to_v100(tmp_path):
    out = tmp_path / "out"
    res = run_cli("run", MAI_GZ, "-o", out, "-t", "2", "--no-plot")
    assert res.exit_code == 0
    for name in ["MSR-1_MAI90kb_mgc_genes.tsv", "MSR-1_MAI90kb_magnetosome_proteins.faa"]:
        assert (out / name).read_bytes() == (EXPECTED / name).read_bytes(), name
    assert (out / "clusters_gbk" / "MSR-1_MAI90kb_MGC.gbk").read_bytes() == \
        (EXPECTED / "MSR-1_MAI90kb_MGC.gbk").read_bytes()


# --------------------------------------------------------------------------
# Fix 1: -s/--score is honoured
# --------------------------------------------------------------------------
MAMK = (
    "MSEGEGQAKNRLFLGIDLGTSHTAVMTSRGKKFLLKSVVGYPKDVIGLKLLGRPYVVGDE"
    "AFEMRSYLDLRYPLQDGVLSEISDRDIEVARHLLTHVVKSAEPGANDEICAVIGVPARAS"
    "GANKALLLKMAQEVVHTALVVSEPFMVGYGLDKLNNTIIVDIGAGTTDICALKGTVPGPE"
    "DQVTLTKAGNYLDERLQNAILERHPELQMNTNVACAVKEQFSFVGARGEAATFEFRAAGK"
    "PVRCDVTESVKIACEALMPDIIESIEILLRSFQPEYQATVLQNIVFAGGGSRIRGLAAYV"
    "KDKLRPFGNADVTCVKDPTFDGCRGALRLAEELPPQYWCQLGDVSGQ"
)


def test_searcher_score_cutoff_semantics():
    gene = PredictedGene("mamK", "c", 1, 3 * len(MAMK), 1, 0, MAMK)
    default = MagnetosomeSearcher()
    assert default.score_cutoff is None
    assert default.min_score_for("MamK") == 25.0          # per-profile trusted_cutoff
    assert "mamK" in default.search_genes([gene])
    strict = MagnetosomeSearcher(score_cutoff=10_000)
    assert strict.min_score_for("MamK") == 10_000
    assert strict.search_genes([gene]) == {}
    assert MagnetosomeSearcher(score_cutoff=5).min_score_for("MamK") == 5.0  # may go below 25


def test_cli_score_changes_hits(tmp_path):
    res_default = run_cli("run", MAI_GZ, "-o", tmp_path / "d", "--no-plot")
    res_high = run_cli("run", MAI_GZ, "-o", tmp_path / "h", "--no-plot", "-s", "300")
    assert res_default.exit_code == res_high.exit_code == 0
    d = tsv_rows(tmp_path / "d" / "MSR-1_MAI90kb_mgc_genes.tsv")
    h = tsv_rows(tmp_path / "h" / "MSR-1_MAI90kb_mgc_genes.tsv")
    assert 0 < len(h) < len(d)
    assert min(float(r.split("\t")[10]) for r in h) >= 300
    log = (tmp_path / "h" / "magcluster.log").read_text()
    assert "score: 300 (explicit" in log


# --------------------------------------------------------------------------
# Fix 2: unique genome IDs and per-genome gene namespaces
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name,expected", [
    ("MSR-1.fna", "MSR-1"),
    ("strainX.v1.fna", "strainX.v1"),
    ("GCF_000513295.1.gbff.gz", "GCF_000513295.1"),
    ("GCA_015228835.1.fna", "GCA_015228835.1"),
    ("genome.FASTA", "genome"),
    ("proteins.faa.gz", "proteins"),
    ("noext", "noext"),
    ("my.strain.data", "my.strain.data"),
])
def test_genome_name_from_path(name, expected):
    assert genome_name_from_path(Path(name)) == expected


def test_same_file_name_in_two_folders(tmp_path, mai_fna):
    for d in ("A", "B"):
        (tmp_path / d).mkdir()
        shutil.copy(mai_fna, tmp_path / d / "genomic.fna")
    out = tmp_path / "out"
    res = run_cli("run", tmp_path / "A" / "genomic.fna", tmp_path / "B" / "genomic.fna", "-o", out, "--no-plot")
    assert res.exit_code == 0
    assert "genomic_2" in flat(res.output) and "already used" in flat(res.output)
    a = tsv_rows(out / "genomic_mgc_genes.tsv")
    b = tsv_rows(out / "genomic_2_mgc_genes.tsv")
    assert len(a) == len(b) == 37
    assert all(r.split("\t")[1] == "genomic" for r in a)
    assert all(r.split("\t")[1] == "genomic_2" for r in b)
    assert sorted((out / "clusters_gbk").iterdir()) == [
        out / "clusters_gbk" / "genomic_2_MGC.gbk", out / "clusters_gbk" / "genomic_MGC.gbk"]


def test_dotted_versions_do_not_overwrite(tmp_path, mai_fna):
    d = tmp_path / "prefix"
    d.mkdir()
    shutil.copy(mai_fna, d / "strainX.v1.fna")
    shutil.copy(mai_fna, d / "strainX.v2.fna")
    out = tmp_path / "out"
    assert run_cli("run", d, "-o", out, "--no-plot").exit_code == 0
    v1 = tsv_rows(out / "strainX.v1_mgc_genes.tsv")
    v2 = tsv_rows(out / "strainX.v2_mgc_genes.tsv")
    assert len(v1) == len(v2) == 37                  # nothing lost to ID collisions
    assert all(r.split("\t")[2].startswith("strainX.v1_") for r in v1)
    assert all(r.split("\t")[2].startswith("strainX.v2_") for r in v2)


def test_identical_gene_ids_in_two_genbank_genomes_stay_separate(tmp_path):
    """Two GenBank files with the same locus tags must not share hits."""
    from magcluster.cli import _assign_genome_names
    names, notices = _assign_genome_names([Path("x/genomic.gbff"), Path("y/genomic.gbff"), Path("z/GENOMIC.gbk")])
    assert names == ["genomic", "genomic_2", "GENOMIC_3"]   # case-insensitive dedupe
    assert len(notices) == 2


# --------------------------------------------------------------------------
# Fix 3: directories with non-genome files; per-genome failures
# --------------------------------------------------------------------------
def test_directory_with_readme_is_processed(tmp_path, mai_fna):
    d = mai_fna.parent
    (d / "README.txt").write_text("notes\n")
    (d / ".hidden.fna").write_text(">x\nACGT\n")
    out = tmp_path / "out"
    res = run_cli("run", d, "-o", out, "--no-plot")
    assert res.exit_code == 0
    assert "Skipped" in flat(res.output) and "README.txt" in flat(res.output)
    assert (out / "MSR-1_MAI90kb_mgc_genes.tsv").exists()


def test_bad_genome_file_does_not_lose_other_results(tmp_path, mai_fna):
    d = mai_fna.parent
    (d / "aaa_broken.fna").write_text("this is not a FASTA file\n")
    (d / "zzz_binary.fna").write_bytes(b"\x00\x01\x02binary")
    out = tmp_path / "out"
    res = run_cli("run", d, "-o", out, "--no-plot")
    assert res.exit_code == 1                         # failures are reported ...
    assert (out / "MSR-1_MAI90kb_mgc_genes.tsv").exists()   # ... but other results are kept
    text = flat(res.output)
    assert "aaa_broken.fna" in text and "zzz_binary.fna" in text
    log = (out / "magcluster.log").read_text()
    assert "Failed input" in log and "aaa_broken.fna" in log


def test_is_genome_file():
    assert is_genome_file(Path("a.fna")) and is_genome_file(Path("a.GBFF.gz")) and is_genome_file(Path("a.pep"))
    assert not is_genome_file(Path("README.txt")) and not is_genome_file(Path("a.gff"))


# --------------------------------------------------------------------------
# Fix 4: spaces / special characters in file names
# --------------------------------------------------------------------------
def test_file_name_with_spaces_and_unicode(tmp_path, mai_fna):
    d = tmp_path / "dir with space" / "中文目录"
    d.mkdir(parents=True)
    src = d / "磁螺菌 MSR-1 (copy).fna"
    shutil.copy(mai_fna, src)
    out = tmp_path / "out dir"
    res = run_cli("run", src, "-o", out, "--no-plot")
    assert res.exit_code == 0
    assert sanitize_genome_name("磁螺菌 MSR-1 (copy)") == "磁螺菌_MSR-1_copy_"
    gbk = out / "clusters_gbk" / "磁螺菌_MSR-1_copy__MGC.gbk"
    recs = list(SeqIO.parse(str(gbk), "genbank"))
    assert recs and all(" " not in r.id for r in recs)
    assert str(src.resolve()) in (out / "magcluster.log").read_text()   # original name traceable


@pytest.mark.skipif(not HAS_CLINKER, reason="clinker not installed")
def test_spaces_html_generated(tmp_path, mai_fna):
    src = tmp_path / "my strain.fna"
    shutil.copy(mai_fna, src)
    res = run_cli("run", src, "-o", tmp_path / "o")
    assert res.exit_code == 0
    assert (tmp_path / "o" / "MGC_align.html").exists()


# --------------------------------------------------------------------------
# Fix 5: dependency floors
# --------------------------------------------------------------------------
def test_dependency_floors():
    text = (REPO / "pyproject.toml").read_text()
    deps = re.search(r"^dependencies = \[(.*?)\]", text, re.S | re.M).group(1)
    assert '"typer>=0.16.0"' in deps
    assert '"pyhmmer>=0.10.1"' in deps
    assert "pandas" not in deps and "pyfamsa" not in deps


# --------------------------------------------------------------------------
# New: run log
# --------------------------------------------------------------------------
def test_run_log_written_and_overwritten(tmp_path):
    out = tmp_path / "out"
    assert run_cli("run", MAI_GZ, "-o", out, "--no-plot", "-th", "4").exit_code == 0
    log = (out / "magcluster.log").read_text()
    from magcluster import __version__
    for needle in [f"MagCluster v{__version__}", "Python:", "pyhmmer", "Command:",
                   "threshold (min distinct non-Feo families): 4", "MSR-1_MAI90kb <-",
                   "MSR-1_MAI90kb_loc1", "Per-genome summary", "exit code 0"]:
        assert needle in log, needle
    assert "[bold" not in log                         # Rich markup stripped
    # --force: log is replaced, not appended
    assert run_cli("run", MAI_GZ, "-o", out, "--no-plot", "--force").exit_code == 0
    log2 = (out / "magcluster.log").read_text()
    assert log2.count("run log") == 1 and "min distinct non-Feo families): 3" in log2


def test_custom_log_path(tmp_path):
    custom = tmp_path / "logs" / "run1.log"
    assert run_cli("run", MAI_GZ, "-o", tmp_path / "out", "--no-plot", "--log", custom).exit_code == 0
    assert custom.exists() and not (tmp_path / "out" / "magcluster.log").exists()


# --------------------------------------------------------------------------
# --force, clinker failures, validation, friendly errors
# --------------------------------------------------------------------------
def test_force_removes_stale_outputs_only(tmp_path, mai_fna):
    out = tmp_path / "out"
    assert run_cli("run", MAI_GZ, "-o", out, "--no-plot", "--static-plot").exit_code == 0
    assert (out / "figures" / "MSR-1_MAI90kb_MGC.svg").exists()
    (out / "my_notes.txt").write_text("keep me")
    # Second run on a genome without any MGC (only 3 kb of the island)
    rec = next(SeqIO.parse(str(mai_fna), "fasta"))
    small = tmp_path / "small.fna"
    small.write_text(f">small\n{str(rec.seq[:3000])}\n")
    assert run_cli("run", small, "-o", out, "--no-plot", "--force").exit_code == 0
    left = sorted(p.name for p in out.iterdir())
    assert "MSR-1_MAI90kb_mgc_genes.tsv" not in left
    assert "clusters_gbk" not in left and "figures" not in left
    assert "my_notes.txt" in left


@pytest.mark.skipif(not HAS_CLINKER, reason="clinker not installed")
def test_clinker_failure_is_an_error(tmp_path):
    res = run_cli("run", MAI_GZ, "-o", tmp_path / "o", "--clinker-opts", "--bogus")
    assert res.exit_code == 1
    text = flat(res.output)
    assert "clinker failed" in text
    assert "Clinker Visualization:" not in text and "✓ Done" not in text
    assert not (tmp_path / "o" / "MGC_align.html").exists()
    assert (tmp_path / "o" / "MSR-1_MAI90kb_mgc_genes.tsv").exists()


@pytest.mark.skipif(not HAS_CLINKER, reason="clinker not installed")
def test_clinker_relative_outputs_go_to_outdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    shutil.copy(MAI_GZ, tmp_path / "a.fna.gz")
    shutil.copy(MAI_GZ, tmp_path / "b.fna.gz")
    res = run_cli("run", "a.fna.gz", "b.fna.gz", "-o", "o", "--clinker-opts", "-mo matrix.tsv")
    assert res.exit_code == 0
    assert (tmp_path / "o" / "matrix.tsv").exists() and not (tmp_path / "matrix.tsv").exists()


@pytest.mark.parametrize("opt", [
    ["-t", "-1"], ["-e", "-1"], ["--max-gap", "-1"], ["-th", "0"], ["-th", "-1"],
    ["-un", "-2"], ["-i", "1.5"], ["-i", "-0.2"], ["--clinker-opts", "-na 'unterminated"],
])
def test_invalid_numeric_and_option_values(tmp_path, opt):
    res = run_cli("run", MAI_GZ, "-o", tmp_path / "o", *opt)
    assert res.exit_code == 2
    assert "Invalid value" in flat(res.output)
    assert not (tmp_path / "o").exists()


def test_threads_above_cpu_count_are_capped(tmp_path):
    res = run_cli("run", MAI_GZ, "-o", tmp_path / "o", "--no-plot", "-t", "100000")
    assert res.exit_code == 0 and "exceeds" in flat(res.output)


@pytest.mark.parametrize("content,needles", [
    (b"", ["no FASTA records"]),
    # Biopython >= 1.85 rejects leading text; older versions return no records
    (b"just some text\nmore\n", ["not a valid FASTA", "no FASTA records"]),
    (b"\x00\x00\x89PNG", ["binary file"]),
])
def test_unreadable_inputs_give_clear_errors(tmp_path, content, needles):
    f = tmp_path / "bad.fna"
    f.write_bytes(content)
    res = run_cli("run", f, "-o", tmp_path / "o")
    assert res.exit_code == 1
    out = flat(res.output)
    assert any(n in out for n in needles)
    # one short line, not Biopython's multi-paragraph explanation
    err = next(line for line in out.split("Error: ")[1:2])
    assert len(err.split("Skipped.")[0]) < 220 and "fasta-pearson" not in out


def test_short_reason_condenses_parser_errors():
    from magcluster.core.parser import _short_reason
    biopython_msg = ("This FASTA file contains comments at the beginning of the file, which are not\n"
                     "allowed by the 'fasta' parser. To parse this file, you have three options: ...")
    assert _short_reason(ValueError(biopython_msg)) == "text found before the first '>' header line"
    long = _short_reason(ValueError("x " * 500))
    assert len(long) <= 160 and long.endswith("...") and "\n" not in long


def test_gff_rejected_clearly(tmp_path):
    f = tmp_path / "genes.gff"
    f.write_text("##gff-version 3\nctg\tpyrodigal\tCDS\t1\t90\t.\t+\t0\tID=1\n")
    with pytest.raises(InputFormatError, match="GFF3 input is not supported"):
        parse_input(f)
    res = run_cli("run", f, "-o", tmp_path / "o")
    assert res.exit_code == 1 and "GFF3 input is not supported" in flat(res.output)


def test_outdir_is_a_file(tmp_path):
    f = tmp_path / "afile"
    f.write_text("x")
    res = run_cli("run", MAI_GZ, "-o", f)
    assert res.exit_code == 1 and "is not a directory" in flat(res.output)


def test_single_mode_on_short_sequence_falls_back(tmp_path, mai_fna):
    rec = next(SeqIO.parse(str(mai_fna), "fasta"))
    short = tmp_path / "short.fna"
    short.write_text(f">s\n{str(rec.seq[40000:56000])}\n")
    res = run_cli("run", short, "-o", tmp_path / "o", "--single", "--no-plot")
    assert res.exit_code == 0
    assert "used metagenome mode instead" in flat(res.output)


def test_duplicate_inputs_are_deduplicated(tmp_path):
    out = tmp_path / "o"
    res = run_cli("run", MAI_GZ, MAI_GZ, "-o", out, "--no-plot")
    assert res.exit_code == 0 and "duplicate input" in flat(res.output)
    assert len(tsv_rows(out / "MSR-1_MAI90kb_mgc_genes.tsv")) == 37


def test_meta_flag_ignored_notice_for_annotated_input(tmp_path):
    faa = tmp_path / "p.faa"
    faa.write_text(">p1 # 1 # 300 # 1 # ID=1;contig=c1\n" + MAMK + "\n")
    res = run_cli("run", faa, "-o", tmp_path / "o", "--meta", "--no-plot")
    assert res.exit_code == 0 and "--meta ignored" in flat(res.output)


# --------------------------------------------------------------------------
# Origin-spanning MGC coordinates
# --------------------------------------------------------------------------
def _hit(gid, ctg, fam, start, idx):
    return MagnetosomeHit(gid, ctg, fam, fam, 100.0, 1e-20, "core", "Alpha", "#E63946",
                          start, start + 899, 1, idx, "M" * 100)


def test_mgc_spanning_origin_on_circular_replicon():
    genes, hits = [], {}
    starts = [1000 + i * 1000 for i in range(5)] + [90_000 + i * 1000 for i in range(5)]
    fams = ["MamA", "MamB", "MamE", None, None, None, None, "MamK", "MamM", "MamO"]
    for i, (st, fam) in enumerate(zip(starts, fams)):
        genes.append(PredictedGene(f"g{i}", "chr", st, st + 899, 1, i, "M" * 100))
        if fam:
            hits[f"g{i}"] = _hit(f"g{i}", "chr", fam, st, i)
    engine = ClusterEngine(flank_genes=0)
    linear = engine.cluster_genome("g", genes, hits, {"chr": 100_000}, {"chr": False})[0]
    assert (linear.start_bp, linear.end_bp, linear.wraps_origin) == (1000, 94899, False)
    circ = engine.cluster_genome("g", genes, hits, {"chr": 100_000}, {"chr": True})[0]
    assert circ.wraps_origin and (circ.start_bp, circ.end_bp) == (92_000, 3899)
    assert circ.span_bp == 100_000 - 92_000 + 1 + 3899
    # membership is unchanged; only the reported span differs
    assert [g.gene_id for g in circ.all_genes] == [g.gene_id for g in linear.all_genes]


def test_cluster_engine_rejects_negative_values():
    with pytest.raises(ValueError):
        ClusterEngine(flank_genes=-1)
    with pytest.raises(ValueError):
        ClusterEngine(min_distinct_families=0)


# --------------------------------------------------------------------------
# Help / start-up
# --------------------------------------------------------------------------
def test_dash_h_alias():
    for args in (["-h"], ["run", "-h"], ["db", "-h"]):
        res = run_cli(*args)
        assert res.exit_code == 0 and re.search(r"--help\s+-h", flat(res.output))


def test_usage_lines_not_blank():
    # typer 0.27 propagates the app-level options_metavar="" to sub-commands,
    # which used to render an empty "Usage:" line for `magcluster db -h`.
    db = flat(run_cli("db", "--help").output)
    assert "Usage: magcluster db [OPTIONS]" in db
    run = flat(run_cli("run", "--help").output)
    assert "Usage: magcluster run [OPTIONS]" in run


def test_clinker_help_shows_default():
    res = run_cli("run", "--clinker-help")
    assert res.exit_code == 0 and "[default: 2]" in flat(res.output)


def test_no_matplotlib_on_import_or_version():
    code = "import sys, magcluster.cli; print('matplotlib' in sys.modules, 'pygenomeviz' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert strip_ansi(out.stdout).strip() == "False False"
    res = subprocess.run([sys.executable, "-c", "from magcluster.cli import app; app()", "--version"],
                         capture_output=True, text=True)
    assert res.returncode == 0
    assert strip_ansi(res.stderr).strip() == ""
    assert strip_ansi(res.stdout).strip().endswith(__import__("magcluster").__version__)
