"""Logging and rich terminal UI helpers for MagCluster."""

import logging
import sys
from rich.console import Console
from rich.logging import RichHandler
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

def setup_logger(level: str = "INFO") -> logging.Logger:
    """Setup and return application logger with rich formatting."""
    log_level = getattr(logging, level.upper(), logging.INFO)
    
    # Configure root logger
    logging.basicConfig(
        level=log_level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[
            RichHandler(
                console=console,
                rich_tracebacks=True,
                show_path=False,
                markup=True
            )
        ]
    )
    
    logger = logging.getLogger("magcluster")
    logger.setLevel(log_level)
    return logger

log = setup_logger()
