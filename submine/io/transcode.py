"""Format detection and transcoding utilities.

Many third-party miners operate on a *native on-disk* format (e.g., gSpan datasets
or SoPaGraMi ``.lg``). For these miners, the most efficient pipeline is:

    user input file  ->  transcode to miner native file  ->  miner runs

This module implements the 'transcode to native' step so the API does not need to
round-trip through an intermediate :class:`~submine.core.graph.Graph` unless the
input itself is not in the miner's native format.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional
import re
from ..core.graph import Graph


class UnknownFormatError(ValueError):
    pass


@dataclass(frozen=True)
class FormatSpec:
    """A lightweight format descriptor."""

    key: str
    suffixes: tuple[str, ...]


# Canonical format keys used across the library.
FMT_LG = "lg"  # SoPaGraMi single-graph format
FMT_GSPAN = "gspan"  # gSpan dataset format (t/v/e ... t#-1)
FMT_GEXF = "gexf"  # NetworkX-readable GEXF
FMT_EDGELIST = "edgelist"  # whitespace-separated u v [label]

_GSPAN_DATA_RE = re.compile(r".*\.data(\.[A-Za-z0-9_-]+)?$")  # matches .data, .data.x, .data.2 etc
_KNOWN_FORMATS: List[FormatSpec] = [
    FormatSpec(FMT_LG, (".lg",)),
    # treat both classic .gspan and Gatech-like *.data / *.data.x as gSpan datasets
    FormatSpec(FMT_GSPAN, (".gspan", ".data", ".data.x")),
    FormatSpec(FMT_GEXF, (".gexf",)),
    FormatSpec(FMT_EDGELIST, (".edgelist", ".txt", ".tsv", ".csv")),
]


def detect_format(path: str | Path) -> str:
    """Detect the most likely input format from a file path.

    Detection is filename-based (not content-based) by design for speed.

    Notes
    -----
    - Gatech-style gSpan datasets frequently use ``*.data`` and ``*.data.<tag>`` where
      ``<tag>`` may be an integer shard index (e.g., ``.data.2``) or an arbitrary token
      (e.g., ``.data.x``). We treat all of these as :data:`FMT_GSPAN`.
    - For ``.txt/.csv/.tsv`` we assume edge-list unless otherwise specified.
    """
    p = Path(path)
    name = p.name.lower()
    suf = p.suffix.lower()

    # Handle multi-suffix gSpan dataset conventions first.
    if _GSPAN_DATA_RE.match(name):
        return FMT_GSPAN

    # Exact suffix matches (single suffix) for known formats.
    for spec in _KNOWN_FORMATS:
        if suf in spec.suffixes:
            return spec.key

    # Fallback for compound suffixes like ".data.x" where Path.suffix == ".x".
    for spec in _KNOWN_FORMATS:
        for ss in spec.suffixes:
            if name.endswith(ss):
                return spec.key

    raise UnknownFormatError(f"Cannot detect graph format from file: {p}")


def load_graphs(path: str | Path, *, fmt: Optional[str] = None) -> List[Graph]:
    """Load graphs from a supported file format."""
    p = Path(path)
    fmt = fmt or detect_format(p)

    if fmt == FMT_GSPAN:
        from .gspan import read_gspan_dataset

        return list(read_gspan_dataset(p))

    if fmt == FMT_LG:
        from .sopagrami import read_lg

        return [read_lg(p)]

    if fmt == FMT_EDGELIST:
        from .common import read_edgelist_dataset

        return list(read_edgelist_dataset(p))

    if fmt == FMT_GEXF:
        from .gexf import read_gexf

        return [read_gexf(p)]

    raise UnknownFormatError(f"Unsupported input format: {fmt}")


def write_graphs(graphs: Iterable[Graph], path: str | Path, *, fmt: str) -> Path:
    """Write graphs to a given native format."""
    p = Path(path)

    if fmt == FMT_GSPAN:
        from .gspan import write_gspan_dataset

        write_gspan_dataset(list(graphs), p)
        return p

    if fmt == FMT_LG:
        from .sopagrami import write_lg

        gs = list(graphs)
        if len(gs) != 1:
            raise ValueError(f".lg expects exactly one graph; got {len(gs)}")
        write_lg(gs[0], p)
        return p

    raise UnknownFormatError(f"Unsupported output format: {fmt}")


def transcode_path(
    src_path: str | Path,
    dst_path: str | Path,
    *,
    dst_fmt: str,
    src_fmt: Optional[str] = None,
) -> Path:
    """Transcode an on-disk graph dataset to another format.

    This function parses the input *once* into in-memory graphs (only when needed)
    and then writes the target native file.
    """
    graphs = load_graphs(src_path, fmt=src_fmt)
    return write_graphs(graphs, dst_path, fmt=dst_fmt)
