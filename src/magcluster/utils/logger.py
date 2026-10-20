"""Logging and rich terminal UI helpers for MagCluster.

Console output is produced with Rich.  Every ``magcluster run`` additionally
records a plain-text run log (``magcluster.log``) through a standard
:class:`logging.FileHandler`: messages shown on the console are mirrored to the
file (with Rich markup removed), and extra detail that is too verbose for the
terminal is written to the file only.
"""

import logging
from pathlib import Path
from typing import Optional
from rich.console import Console
from rich.logging import RichHandler
from rich.text import Text
from rich.theme import Theme

custom_theme = Theme({
    "info": "cyan",
    "warning": "yellow",
    "error": "bold red",
    "success": "bold green",
    "highlight": "bold magenta",
})

console = Console(theme=custom_theme)
err_console = Console(stderr=True, theme=custom_theme)

LOGGER_NAME = "magcluster"
#: Logger that only writes to the run-log file (never to the console).
RUNLOG_NAME = "magcluster.runlog"


def _strip_markup(msg: str) -> str:
    try:
        return Text.from_markup(msg).plain
    except Exception:
        return msg


class PlainFormatter(logging.Formatter):
    """File formatter: timestamp + level + message with Rich markup removed."""

    def __init__(self):
        super().__init__(fmt="%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

    def format(self, record: logging.LogRecord) -> str:
        record = logging.makeLogRecord(record.__dict__)
        record.msg = _strip_markup(record.getMessage())
        record.args = None
        return super().format(record)


def setup_logger(level: str = "INFO") -> logging.Logger:
    """Return the application logger with Rich console formatting.

    Only the ``magcluster`` logger is configured (the root logger is left alone),
    so INFO messages of third-party libraries such as matplotlib never reach the
    terminal.
    """
    log_level = getattr(logging, level.upper(), logging.INFO)
    logger = logging.getLogger(LOGGER_NAME)
    if not any(isinstance(h, RichHandler) for h in logger.handlers):
        handler = RichHandler(
            console=console,
            rich_tracebacks=True,
            show_path=False,
            markup=True,
        )
        handler.setFormatter(logging.Formatter("%(message)s", datefmt="[%X]"))
        logger.addHandler(handler)
    logger.setLevel(log_level)
    logger.propagate = False

    runlog = logging.getLogger(RUNLOG_NAME)
    runlog.setLevel(logging.DEBUG)
    runlog.propagate = False
    return logger


log = setup_logger()
runlog = logging.getLogger(RUNLOG_NAME)


def attach_log_file(path: Path) -> logging.FileHandler:
    """Start writing the run log to *path* (overwritten).  Returns the handler."""
    path.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    handler.setFormatter(PlainFormatter())
    handler.setLevel(logging.DEBUG)
    logging.getLogger(LOGGER_NAME).addHandler(handler)
    runlog.addHandler(handler)
    return handler


def detach_log_file(handler: Optional[logging.FileHandler]):
    """Stop writing to the run-log file and close it."""
    if handler is None:
        return
    for name in (LOGGER_NAME, RUNLOG_NAME):
        logging.getLogger(name).removeHandler(handler)
    handler.close()


def say(msg: str = "", level: int = logging.INFO):
    """Print *msg* to the console (Rich markup allowed) and mirror it to the run log."""
    console.print(msg)
    if msg.strip():
        runlog.log(level, msg.strip())


def detail(msg: str, level: int = logging.INFO):
    """Write *msg* to the run log only (not shown on the console)."""
    runlog.log(level, msg)


def quiet_matplotlib():
    """Silence matplotlib's font-cache / font_manager chatter before importing it."""
    for name in ("matplotlib", "matplotlib.font_manager", "fontTools"):
        logging.getLogger(name).setLevel(logging.ERROR)
