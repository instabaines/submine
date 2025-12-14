
from __future__ import annotations

from pathlib import Path
from typing import Dict, Hashable

from ..core.graph import Graph


def write_lg(graph: Graph, path: str | Path, directed: bool = False) -> None:
    """
    Write a single Graph to SoPaGraMi's .lg format.

    - Nodes are reindexed to 0..n-1 internally.
    - Node labels come from graph.node_labels (fallback to string of node id).
    - Edge labels come from graph.edge_labels (fallback to empty string).
    """
    path = Path(path)

    # Map original node ids -> contiguous [0..n-1]
    node_ids = list(graph.nodes)
    id_map: Dict[Hashable, int] = {nid: i for i, nid in enumerate(node_ids)}

    node_labels = graph.node_labels or {}
    edge_labels = graph.edge_labels or {}

    with path.open("w") as f:
        # vertices
        for nid in node_ids:
            idx = id_map[nid]
            lbl = node_labels.get(nid, str(nid))
            f.write(f"v {idx} {lbl}\n")

        # edges
        for (u_orig, v_orig) in graph.edges:
            u = id_map[u_orig]
            v = id_map[v_orig]

            # SoPaGraMi supports directed; if undirected we still write one edge
            lbl = edge_labels.get((u_orig, v_orig)) or edge_labels.get((v_orig, u_orig)) or ""
            if lbl == "":
                f.write(f"e {u} {v}\n")
            else:
                f.write(f"e {u} {v} {lbl}\n")
