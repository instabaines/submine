"""Common I/O helpers used by algorithm wrappers.

Functions defined here are not tied to a specific algorithm. They
perform tasks such as ensuring that directories exist or creating
temporary working directories.
"""

from __future__ import annotations

import os
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Iterator

from pathlib import Path
from typing import List, Optional, Tuple, Any

from ..core.graph import Graph
from ..errors import SubmineInputError
from ..utils.checks import iter_text_lines


def ensure_dir(path: str | Path) -> Path:
    """Ensure that the directory at ``path`` exists.

    Parameters
    ----------
    path: str or pathlib.Path
        Directory path to create if it does not exist.

    Returns
    -------
    pathlib.Path
        The Path instance corresponding to ``path``.
    """
    p = Path(path)
    if not p.exists():
        p.mkdir(parents=True, exist_ok=True)
    return p


@contextmanager
def temporary_directory() -> Iterator[Path]:
    """Context manager yielding a temporary directory as a Path.

    The directory and its contents are removed on exit. Useful when
    writing temporary files for external algorithms.
    """
    with TemporaryDirectory() as tmp:
        yield Path(tmp)


def _maybe_int(x: str) -> Any:
    """Try to parse as int, fallback to string."""
    try:
        return int(x)
    except ValueError:
        return x


def read_edgelist_dataset(path: str | Path) -> List[Graph]:
    """
    Read an edge-list dataset file and return a list of Graph objects.

    Supported formats
    -----------------
    1) Single graph (no 't # gid' headers):

        u v
        u v label

       - Nodes are inferred from all endpoints.
       - Edge labels are used if a third column is present on ANY line.

    2) Multiple graphs (with gSpan-like headers):

        t # 0
        u v
        u v label
        t # 1
        u v
        ...

       - Each 't # gid' starts a new graph.
       - Nodes are inferred per graph.

    Returns
    -------
    List[Graph]
        A list of Graph objects constructed from the file.
    """
    path = Path(path)

    graphs: List[Graph] = []

    current_edges: Optional[List[Tuple[Any, Any]]] = None
    current_edge_labels: Optional[dict[Tuple[Any, Any], Any]] = None
    any_labels_in_current = False

    def flush_current_graph():
        nonlocal current_edges, current_edge_labels, any_labels_in_current
        if current_edges is None:
            return

        # Build node set from edges
        nodes_set = set()
        for u, v in current_edges:
            nodes_set.add(u)
            nodes_set.add(v)
        # Deterministic node ordering for reproducible transcoding/writing.
        # Edge lists may mix ints/strings; sort by type name then by string value.
        nodes = sorted(nodes_set, key=lambda x: (type(x).__name__, str(x)))

        edge_labels = current_edge_labels if any_labels_in_current else None

        graphs.append(
            Graph(
                nodes=nodes,
                edges=current_edges,
                node_labels=None,        # edge-list doesn't carry node labels
                edge_labels=edge_labels,
            )
        )
        current_edges = None
        current_edge_labels = None
        any_labels_in_current = False

    for raw_line in iter_text_lines(path):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        # Accept both whitespace-separated and comma-separated edge lists.
        if "," in line:
            line = line.replace(",", " ")
        parts = line.split()
        rec_type = parts[0]

        # Multi-graph header: t # gid
        if rec_type == "t" and len(parts) >= 3 and parts[1] == "#":
            # flush previous graph, start new
            flush_current_graph()
            current_edges = []
            current_edge_labels = {}
            any_labels_in_current = False
            # we ignore gid value; it's just a marker
            continue

        # Otherwise, treat as an edge line: u v [label]
        if current_edges is None:
            # No graph header seen yet: assume single-graph file
            current_edges = []
            current_edge_labels = {}
            any_labels_in_current = False

        if len(parts) < 2:
            raise ValueError(f"Malformed edge line: {line!r}")

        u = _maybe_int(parts[0])
        v = _maybe_int(parts[1])

        if len(parts) >= 3:
            lbl = _maybe_int(parts[2])
            any_labels_in_current = True
            current_edges.append((u, v))
            current_edge_labels[(u, v)] = lbl  # type: ignore[index]
        else:
            current_edges.append((u, v))

    # Flush final graph (single-graph files without a trailing header)
    flush_current_graph()
    return graphs