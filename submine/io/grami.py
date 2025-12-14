"""Utilities for GraMi algorithm input/output.

GraMi accepts graphs in a simple text format where each graph is
represented by a sequence of lines. There is no official Python
interface, so this module provides helper functions to emit a
compatible format. The current implementation is rudimentary and may
need adjustments depending on the GraMi version you are using.
"""

from __future__ import annotations

from typing import List

from ..core.graph import Graph


def write_grami_input(graphs: List[Graph], file_path: str) -> None:
    """Write graphs to disk in a basic GraMi input format.

    Each graph is written as a series of lines with the following
    structure::

        t # <graph_id>
        v <node_id> <label>
        e <src> <dst> <label>

    Similar to gSpan, node identifiers are renumbered per graph. This
    format is accepted by some versions of GraMi, but you may need to
    consult the documentation for your particular release.

    Parameters
    ----------
    graphs: List[Graph]
        Graphs to serialise.
    file_path: str
        Path to the output file. Parent directories should already
        exist.
    """
    with open(file_path, "w", encoding="utf-8") as f:
        for gid, g in enumerate(graphs):
            f.write(f"t # {gid}\n")
            id_map = {node.id: idx for idx, node in enumerate(g.nodes())}
            for node in g.nodes():
                nid = id_map[node.id]
                label = node.label or str(node.id)
                f.write(f"v {nid} {label}\n")
            for u, v, edge_label in g.iter_edges():
                u_idx = id_map[u]
                v_idx = id_map[v]
                label = edge_label if edge_label is not None else "1"
                f.write(f"e {u_idx} {v_idx} {label}\n")