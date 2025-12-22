"""Environment and safety checks for submine.

This module centralizes lightweight validation utilities that protect the
library from common failure modes (corrupt inputs, missing binaries) and
from avoidable abuse (e.g., attempting to load arbitrarily large files into
memory, unsafe subprocess execution defaults).

These are not meant to be a sandbox. They are meant to provide sensible,
defensive defaults for a publishable OSS library.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Optional


DEFAULT_MAX_INPUT_MB = int(os.getenv("SUBMINE_MAX_INPUT_MB", "128"))
DEFAULT_MAX_INPUT_BYTES = DEFAULT_MAX_INPUT_MB * 1024 * 1024


DEFAULT_MAX_LINES = int(os.getenv("SUBMINE_MAX_LINES", "5000000"))  # 5M lines
DEFAULT_MAX_LINE_BYTES = int(os.getenv("SUBMINE_MAX_LINE_BYTES", "1048576"))  # 1 MiB per line


def iter_text_lines(
    path: str | Path,
    *,
    encoding: str = "utf-8",
    max_lines: int | None = None,
    max_line_bytes: int | None = None,
):
    """Yield decoded lines from *path* with hard limits.

    This is intended for parsers that stream large graph files. It protects
    against:
      - extremely long lines (often accidental corruption or malicious inputs)
      - unbounded files that could consume excessive CPU time

    Notes:
      - Lines are decoded with ``errors='replace'`` to avoid UnicodeDecodeError.
      - The returned lines are stripped of trailing ``\n``.
    """
    # Resolve limits at call time so test suites (and embedding apps) can
    # override them via environment variables without requiring a reload.
    if max_lines is None:
        max_lines = int(os.getenv("SUBMINE_MAX_LINES", str(DEFAULT_MAX_LINES)))
    if max_line_bytes is None:
        max_line_bytes = int(os.getenv("SUBMINE_MAX_LINE_BYTES", str(DEFAULT_MAX_LINE_BYTES)))

    p = assert_regular_file(path)
    count = 0
    with p.open("rb") as f:
        for raw in f:
            count += 1
            if count > max_lines:
                from ..errors import ResourceLimitError
                raise ResourceLimitError(
                    f"Refusing to parse {p}: exceeds max line count limit ({max_lines}). "
                    "Set SUBMINE_MAX_LINES to increase the limit if you trust this input."
                )
            if len(raw) > max_line_bytes:
                from ..errors import ResourceLimitError
                raise ResourceLimitError(
                    f"Refusing to parse {p}: line {count} exceeds max line length ({max_line_bytes} bytes). "
                    "Set SUBMINE_MAX_LINE_BYTES to increase the limit if you trust this input."
                )
            yield raw.decode(encoding, errors="replace").rstrip("\n")

def is_tool_available(name: str) -> bool:
    """Return True if a given executable exists on the system PATH."""
    return shutil.which(name) is not None


def assert_regular_file(path: str | Path, *, must_exist: bool = True) -> Path:
    """Validate that *path* points to a regular file.

    We reject directories and (by default) require existence. We also resolve
    symlinks to avoid surprises.
    """
    p = Path(path).expanduser()
    if must_exist and not p.exists():
        raise FileNotFoundError(f"File does not exist: {p}")
    if must_exist and not p.is_file():
        raise ValueError(f"Expected a regular file path, got: {p}")
    # Resolve to eliminate '..' segments and follow symlinks.
    return p.resolve()


def assert_file_size_under(path: str | Path, *, max_bytes: int = DEFAULT_MAX_INPUT_BYTES) -> None:
    """Raise if the file exceeds *max_bytes*.

    This is primarily used to protect code paths that necessarily read the full
    file into memory (e.g., certain bindings).
    """
    p = Path(path)
    try:
        size = p.stat().st_size
    except OSError as e:
        raise OSError(f"Unable to stat file: {p}") from e
    if size > max_bytes:
        from ..errors import ResourceLimitError
        raise ResourceLimitError(
            f"Refusing to load {p} ({size} bytes): exceeds configured limit of {max_bytes} bytes. "
            "Set SUBMINE_MAX_INPUT_MB to increase the limit if you trust this input."
        )


def safe_read_text(path: str | Path, *, encoding: str = "utf-8", max_bytes: int = DEFAULT_MAX_INPUT_BYTES) -> str:
    """Read a text file with a hard cap on bytes."""
    p = assert_regular_file(path)
    assert_file_size_under(p, max_bytes=max_bytes)
    return p.read_text(encoding=encoding, errors="replace")
