"""Input/Output utilities for submine.

This package contains helper functions to serialise graphs to the input
formats expected by different subgraph mining algorithms and to parse
their outputs back into :class:`~submine.core.graph.Graph` objects.
"""

from .common import ensure_dir  # noqa: F401
from .gspan import write_gspan_dataset  # noqa: F401
from .transcode import (  # noqa: F401
    detect_format,
    load_graphs,
    transcode_path,
    FMT_EDGELIST,
    FMT_GEXF,
    FMT_GSPAN,
    FMT_LG,
)

__all__ = [
    "ensure_dir",
    "write_gspan_dataset",
    "detect_format",
    "load_graphs",
    "transcode_path",
    "FMT_EDGELIST",
    "FMT_GEXF",
    "FMT_GSPAN",
    "FMT_LG",
]