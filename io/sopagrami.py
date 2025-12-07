"""I/O utilities for SoPaGraMi algorithm.

The SoPaGraMi algorithm may require input in its own specific format
distinct from gSpan or GraMi. Since there is no universally agreed
format, the functions provided here are placeholders illustrating
where such logic would live. Developers integrating a concrete
implementation should modify these functions accordingly.
"""

from __future__ import annotations

from typing import List

from ..core.graph import Graph


def write_sopagrami_input(graphs: List[Graph], file_path: str) -> None:
    """Write graphs for SoPaGraMi input (placeholder implementation)."""
    with open(file_path, "w", encoding="utf-8") as f:
        for gid, g in enumerate(graphs):
            f.write(f"# graph {gid}\n")
            # Write vertices: id and label
            for node in g.nodes():
                label = node.label or str(node.id)
                f.write(f"v {node.id} {label}\n")
            # Write edges with optional labels
            for u, v, edge_label in g.iter_edges():
                label = edge_label if edge_label is not None else "1"
                f.write(f"e {u} {v} {label}\n")