"""Logging utilities for submine.

This module configures a simple hierarchical logging setup that
suppresses log output by default. Users can enable verbose logging
per algorithm by passing ``verbose=True`` to the constructor. All
modules within submine should obtain their logger via
:func:`get_logger` rather than calling :func:`logging.getLogger`
directly.
"""

from __future__ import annotations

import logging
from typing import Optional


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """Return a logger with a default configuration.

    If called for the first time this function sets up a root logger
    with a basic configuration that logs messages with level WARNING and
    above to stderr. Subsequent calls return child loggers that
    propagate messages to the root.

    Parameters
    ----------
    name: str, optional
        Name of the logger. If omitted, a root logger is returned.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    if not logging.getLogger().handlers:
        # Configure root logger once
        logging.basicConfig(
            level=logging.WARNING,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
    return logging.getLogger(name)