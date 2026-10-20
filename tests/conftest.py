"""Shared test set-up.

The CLI prints through Rich consoles whose width and colour depend on the
environment (``COLUMNS``, ``FORCE_COLOR``, ``GITHUB_ACTIONS``...).  Tests check
*content*, so MagCluster's own consoles are given a fixed, wide width here and
assertions strip ANSI escape codes (see ``strip_ansi`` in the test modules).
"""

import pytest

from magcluster.utils import logger as _mc_logger


@pytest.fixture(autouse=True)
def _wide_magcluster_console():
    consoles = (_mc_logger.console, _mc_logger.err_console)
    saved = [c._width for c in consoles]
    for c in consoles:
        c.width = 200
    yield
    for c, w in zip(consoles, saved):
        c._width = w
