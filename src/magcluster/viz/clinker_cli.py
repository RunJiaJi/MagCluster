"""Standalone ``magcluster clinker`` subcommand helpers and execution."""

import gzip
import logging
import platform
import shlex
import shutil
import sys
import tempfile
import time
import traceback
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import typer
from rich.markup import escape
from .. import __version__
from ..core.parser import (
    GENBANK_EXTENSIONS,
    InputFormatError,
    _strip_gz,
    known_extension,
    sanitize_genome_name,
)
from ..utils.logger import attach_log_file, console, detach_log_file, detail, say
from .clinker_export import (
    ClinkerError,
    extract_hits_from_genbank,
    parse_clinker_opts,
    run_clinker_batches,
)

MAX_CONSOLE_HTML_LIST = 10


def is_genbank_path(path: Path) -> bool:
    """True if *path* has a recognised GenBank extension (optionally ``.gz``)."""
    return known_extension(path) in GENBANK_EXTENSIONS


def gbk_cluster_stem(path: Path) -> str:
    """Cluster stem from a GenBank path (strips ``.gz`` and GenBank extension)."""
    name = _strip_gz(path.name)
    ext = known_extension(path)
    if ext and name.lower().endswith(ext):
        name = name[: -len(ext)]
    return sanitize_genome_name(name or path.stem)


def collect_clinker_inputs(
    inputs: Optional[List[Path]],
    list_file: Optional[Path] = None,
) -> Tuple[List[Path], List[Tuple[Path, str]], List[Path], List[str], List[Tuple[str, str]]]:
    """Collect and order GenBank files from CLI paths and/or ``--list``.

    Returns ``(files, skipped, duplicates, notices, missing_from_list)``.
    """
    candidates: List[Path] = []
    skipped: List[Tuple[Path, str]] = []
    notices: List[str] = []

    for inp in inputs or []:
        if inp.is_dir():
            direct_gbks: List[Path] = []
            dir_other: List[Tuple[Path, str]] = []
            for p in sorted(inp.iterdir(), key=lambda x: x.name):
                if p.name.startswith(".") or not p.is_file():
                    continue
                if is_genbank_path(p):
                    direct_gbks.append(p)
                else:
                    dir_other.append((p, "not a GenBank (.gbk/.gbff) file"))

            sub_gbk_dir = inp / "clusters_gbk"
            if not direct_gbks and sub_gbk_dir.is_dir():
                sub_gbks = [
                    p for p in sorted(sub_gbk_dir.iterdir(), key=lambda x: x.name)
                    if not p.name.startswith(".") and p.is_file() and is_genbank_path(p)
                ]
                if sub_gbks:
                    notices.append(f"Auto-detected sliced GenBank directory: {sub_gbk_dir}")
                    candidates.extend(sub_gbks)
                    continue

            candidates.extend(direct_gbks)
            skipped.extend(dir_other)
        else:
            candidates.append(inp)

    missing_entries: List[Tuple[str, str]] = []
    if list_file is not None:
        # Build lookup table over discovered directory/file candidates (if any)
        by_key: Dict[str, Path] = {}
        for p in candidates:
            stem = gbk_cluster_stem(p)
            for k in (str(p), str(p.resolve()), p.name, _strip_gz(p.name), stem):
                by_key.setdefault(k.lower(), p)
            if stem.lower().endswith("_mgc"):
                by_key.setdefault(stem[:-4].lower(), p)

        ordered: List[Path] = []
        raw_lines = list_file.read_text(encoding="utf-8", errors="ignore").splitlines()
        for raw in raw_lines:
            entry = raw.strip()
            if not entry or entry.startswith("#"):
                continue
            matched = by_key.get(entry.lower())
            if matched is not None:
                ordered.append(matched)
                continue
            cand_cwd = Path(entry)
            cand_rel = list_file.parent / entry
            if cand_cwd.is_file():
                ordered.append(cand_cwd)
            elif cand_rel.is_file():
                ordered.append(cand_rel)
            else:
                missing_entries.append((entry, f"listed entry '{entry}' not found"))
        candidates = ordered

    files: List[Path] = []
    duplicates: List[Path] = []
    seen = set()
    for p in candidates:
        key = p.resolve()
        if key in seen:
            duplicates.append(p)
            continue
        seen.add(key)
        files.append(p)
    return files, skipped, duplicates, notices, missing_entries


def _existing_html_outputs(work_dir: Path, base_html: Path) -> List[Path]:
    """Find existing target HTML or batched HTML files matching *base_html*."""
    found: List[Path] = []
    if base_html.exists() and (base_html.is_file() or base_html.is_symlink()):
        found.append(base_html)
    if work_dir.is_dir():
        suffix = base_html.suffix or ".html"
        stem = base_html.stem if base_html.suffix else base_html.name
        for p in sorted(work_dir.glob(f"{stem}_batch*{suffix}")):
            if (p.is_file() or p.is_symlink()) and p not in found:
                found.append(p)
    return found


def run_clinker_subcommand(
    *,
    inputs: Optional[List[Path]],
    list_file: Optional[Path],
    output: Path,
    batch_size: int,
    identity: float,
    threads: int,
    clinker_opts: Optional[str],
    force: bool,
    log_path: Optional[Path],
    available_cpus: int,
    dep_versions: Dict[str, str],
) -> None:
    """Validate inputs, prepare GenBank slices, and run standalone Clinker plotting."""
    start_time = time.time()

    try:
        parse_clinker_opts(clinker_opts)
    except ValueError as e:
        raise typer.BadParameter(str(e), param_hint="'--clinker-opts'")

    if not inputs and list_file is None:
        raise typer.BadParameter(
            "Provide at least one GenBank file/directory or --list FILE.",
            param_hint="'INPUTS'",
        )

    is_html_target = output.suffix.lower() in {".html", ".htm"}
    if is_html_target:
        if output.exists() and output.is_dir():
            console.print(
                f"[bold red]Error:[/bold red] Output path '[bold]{escape(str(output))}[/bold]' is an existing directory, "
                f"not an HTML file."
            )
            raise typer.Exit(code=1)
        work_dir = output.parent if str(output.parent) != "" else Path(".")
        base_html = output
    else:
        if output.exists() and not output.is_dir():
            console.print(
                f"[bold red]Error:[/bold red] Output path '[bold]{escape(str(output))}[/bold]' exists and is not a directory."
            )
            raise typer.Exit(code=1)
        work_dir = output
        base_html = work_dir / "MGC_align.html"

    console.print(f"[bold cyan]=== MagCluster v{__version__} Clinker Visualization ===[/bold cyan]")
    console.print(f"[dim]Output target: {escape(str((base_html if is_html_target else work_dir).resolve()))}[/dim]\n")

    expanded_files, skipped_files, duplicate_files, notices, missing_entries = collect_clinker_inputs(
        inputs, list_file=list_file
    )

    if not expanded_files and not missing_entries:
        for p, reason in skipped_files:
            console.print(f"[yellow]! Skipped {escape(str(p))}: {reason}.[/yellow]")
        console.print("[bold red]No GenBank (.gbk/.gbff) files found![/bold red]")
        raise typer.Exit(code=1)

    existing_htmls = _existing_html_outputs(work_dir, base_html)
    removed: List[Path] = []
    if existing_htmls:
        if not force:
            shown = ", ".join(escape(str(p)) for p in existing_htmls[:3])
            console.print(
                f"[bold red]Error:[/bold red] Clinker HTML output already exists ({shown}).\n"
                f"Please choose another output path or pass [bold cyan]--force[/bold cyan] (or [bold cyan]-f[/bold cyan]) to overwrite."
            )
            raise typer.Exit(code=1)
        for p in existing_htmls:
            p.unlink()
            removed.append(p)

    try:
        work_dir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        console.print(f"[bold red]Error:[/bold red] cannot create directory '{escape(str(work_dir))}': {e.strerror or e}")
        raise typer.Exit(code=1)

    log_file = log_path if log_path is not None else work_dir / "magcluster_clinker.log"
    try:
        log_handler = attach_log_file(log_file)
    except OSError as e:
        console.print(f"[bold red]Error:[/bold red] cannot write log file '{escape(str(log_file))}': {e.strerror or e}")
        raise typer.Exit(code=1)

    try:
        eff_threads = available_cpus if threads == 0 else min(threads, available_cpus)
        detail(f"MagCluster v{__version__} clinker log")
        detail(f"Started: {datetime.now().astimezone().isoformat(timespec='seconds')}")
        detail(f"Command: {Path(sys.argv[0]).name} {shlex.join(sys.argv[1:])}")
        detail(f"Working directory: {Path.cwd()}")
        detail(f"Python: {platform.python_version()} ({sys.executable}) on {platform.platform()}")
        detail("Dependencies: " + ", ".join(f"{k} {v}" for k, v in dep_versions.items()))
        detail("Effective parameters:")
        for k, v in [
            ("inputs", ", ".join(str(p) for p in (inputs or [])) or "-"),
            ("list_file", str(list_file.resolve()) if list_file else "-"),
            ("output", str((base_html if is_html_target else work_dir).resolve())),
            ("log", str(Path(log_file).resolve())),
            ("threads", f"{eff_threads} (0 = all {available_cpus} available cores)" if threads == 0 else str(eff_threads)),
            ("batch_size", f"{batch_size} (0 = single HTML)" if batch_size == 0 else str(batch_size)),
            ("identity", f"{identity:g}"),
            ("clinker_opts", clinker_opts or "-"),
            ("force", str(force)),
        ]:
            detail(f"  {k}: {v}")
        if removed:
            detail("--force: removed existing HTML file(s): " + ", ".join(str(p) for p in removed))

        for msg in notices:
            say(f"[dim]Note: {escape(msg)}[/dim]")
        for p, reason in skipped_files:
            say(f"[yellow]! Skipped {escape(str(p))}: {reason}.[/yellow]", logging.WARNING)
        for p in duplicate_files:
            say(f"[yellow]! Skipped duplicate input {escape(str(p))} (already listed).[/yellow]", logging.WARNING)
        if threads > available_cpus:
            say(f"[yellow]! --threads {threads} exceeds the {available_cpus} available CPU cores; using {available_cpus}.[/yellow]",
                logging.WARNING)

        failed: List[Tuple[str, str]] = []
        for entry, reason in missing_entries:
            failed.append((entry, reason))
            say(f"  [bold red]Error:[/bold red] {escape(reason)}.", logging.ERROR)

        if not expanded_files:
            say(f"\n[bold red]✗ No valid GenBank inputs to plot ({len(failed)} failed).[/bold red] See the log: {escape(str(log_file))}",
                logging.ERROR)
            raise typer.Exit(code=1)

        say(f"Discovered [bold]{len(expanded_files)}[/bold] GenBank cluster file(s).")

        hits_by_cluster: Dict[str, Dict] = {}
        staged_gbks: List[Path] = []
        used_stems: Dict[str, Path] = {}

        with tempfile.TemporaryDirectory(prefix="magcluster_clinker_") as tmp_dir_str:
            tmp_dir = Path(tmp_dir_str)
            for p in expanded_files:
                try:
                    c_hits = extract_hits_from_genbank(p, cpus=eff_threads)
                except InputFormatError as e:
                    failed.append((str(p), str(e)))
                    say(f"  [bold red]Error:[/bold red] cannot read {escape(p.name)}: {escape(str(e))}. Skipped.",
                        logging.ERROR)
                    continue
                except Exception as e:
                    failed.append((str(p), f"{type(e).__name__}: {e}"))
                    say(f"  [bold red]Error:[/bold red] failed to read {escape(p.name)} "
                        f"({escape(type(e).__name__)}: {escape(str(e))}). Skipped.",
                        logging.ERROR)
                    detail(traceback.format_exc(), logging.ERROR)
                    continue

                raw_stem = gbk_cluster_stem(p)
                unique_stem = raw_stem
                n = 2
                while unique_stem.lower() in used_stems:
                    unique_stem = f"{raw_stem}_{n}"
                    n += 1
                if unique_stem != raw_stem:
                    say(
                        f"[yellow]! Cluster name '{escape(raw_stem)}' is already used by "
                        f"{escape(str(used_stems[raw_stem.lower()]))}; {escape(str(p))} will be plotted as "
                        f"'{escape(unique_stem)}'.[/yellow]",
                        logging.WARNING,
                    )
                used_stems[unique_stem.lower()] = p
                detail(f"  {unique_stem} <- {p.resolve()} ({len(c_hits)} magnetosome gene(s))")

                is_gz = str(p).lower().endswith(".gz")
                needs_staging = is_gz or (unique_stem != p.stem)
                if needs_staging:
                    staged_path = tmp_dir / f"{unique_stem}.gbk"
                    if is_gz:
                        with gzip.open(p, "rt", encoding="utf-8", errors="ignore") as src_f:
                            staged_path.write_text(src_f.read(), encoding="utf-8")
                    else:
                        shutil.copyfile(p, staged_path)
                    staged_gbks.append(staged_path)
                else:
                    staged_gbks.append(p)
                hits_by_cluster[unique_stem] = c_hits

            if not staged_gbks:
                say(f"\n[bold red]✗ No valid GenBank files could be plotted ({len(failed)} failed).[/bold red] "
                    f"See the log: {escape(str(log_file))}", logging.ERROR)
                raise typer.Exit(code=1)

            clinker_failed = False
            clinker_htmls: List[Path] = []

            def _announce(b_idx: int, n_batches: int, target_html: Path, n_cls: int):
                if n_batches == 1:
                    say(f"[dim]-> Running Clinker on {n_cls} cluster(s): {target_html.name}...[/dim]")
                else:
                    say(f"[dim]-> Running Clinker [batch {b_idx}/{n_batches}, {n_cls} clusters]: {target_html.name}...[/dim]")

            try:
                clinker_htmls = run_clinker_batches(
                    gbk_files=staged_gbks,
                    base_html=base_html,
                    hits={},
                    batch_size=batch_size,
                    jobs=eff_threads,
                    identity=identity,
                    clinker_opts=clinker_opts,
                    hits_by_cluster=hits_by_cluster,
                    cwd=work_dir,
                    raise_on_error=True,
                    on_batch_start=_announce,
                )
            except ClinkerError as e:
                clinker_failed = True
                say(f"  [bold red]Error:[/bold red] {escape(str(e))}", logging.ERROR)

        elapsed = time.time() - start_time
        if failed or clinker_failed:
            problems = []
            if failed:
                problems.append(f"{len(failed)} input(s) failed")
            if clinker_failed:
                problems.append("Clinker failed")
            say(f"\n[bold yellow]! Finished with errors in {elapsed:.1f}s ({'; '.join(problems)}).[/bold yellow]",
                logging.WARNING)
        else:
            say(f"\n[bold green]✓ Done in {elapsed:.1f}s![/bold green]")

        if clinker_htmls and not clinker_failed:
            if len(clinker_htmls) <= MAX_CONSOLE_HTML_LIST:
                for hp in clinker_htmls:
                    say(f"  • Clinker Visualization: [underline]{escape(str(hp))}[/underline]")
            else:
                say(f"  • Clinker Visualization: {len(clinker_htmls)} HTML files "
                    f"({escape(clinker_htmls[0].name)} .. {escape(clinker_htmls[-1].name)}) in [underline]{escape(str(work_dir))}[/underline]")
                for hp in clinker_htmls:
                    detail(f"  • Clinker Visualization: {hp}")
        console.print(f"  • Run Log: [underline]{escape(str(log_file))}[/underline]")

        for item, reason in failed:
            detail(f"Failed input: {item}: {reason}", logging.ERROR)
        code = 1 if (failed or clinker_failed) else 0
        detail(f"Finished: {datetime.now().astimezone().isoformat(timespec='seconds')} "
               f"(elapsed {elapsed:.1f}s, exit code {code})")
        if code:
            raise typer.Exit(code=code)
    except BaseException as e:
        if not isinstance(e, (typer.Exit, SystemExit)):
            detail("Clinker command aborted:\n" + traceback.format_exc(), logging.ERROR)
        raise
    finally:
        detach_log_file(log_handler)
