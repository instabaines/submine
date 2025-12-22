
from __future__ import annotations

from pathlib import Path
from typing import Dict, Hashable, Optional, Tuple

from ..core.graph import Graph
from ..errors import SubmineInputError
from ..utils.checks import iter_text_lines


def read_lg(path: str | Path) -> Graph:
    """Read a single SoPaGraMi/gSpan-style ``.lg`` file into :class:`~submine.core.graph.Graph`.

    Supports edge lines of the form:

    - ``e u v``
    - ``e u v label``
    - ``e u v label weight``
    - ``e u v weight`` (rare; treated as unlabeled edge with weight)

    The reader is streaming and suitable for large graphs.
    """
    path = Path(path)

    nodes: list[int] = []
    node_labels: dict[int, str] = {}
    edges: list[Tuple[int, int]] = []
    edge_labels: dict[Tuple[int, int], str] = {}
    edge_weights: dict[Tuple[int, int], float] = {}

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            rec = parts[0]

            if rec == "t":
                # dataset marker (ignored); SoPaGraMi uses a single graph in practice
                continue

            if rec == "v":
                if len(parts) < 3:
                    raise SubmineInputError(f"Malformed vertex line: {line!r}")
                vid = int(parts[1])
                lbl = parts[2]
                if vid not in node_labels:
                    nodes.append(vid)
                node_labels[vid] = lbl
                continue

            if rec == "e":
                if len(parts) < 3:
                    raise SubmineInputError(f"Malformed edge line: {line!r}")
                u = int(parts[1])
                v = int(parts[2])
                edges.append((u, v))

                # Parse optional label / weight. We accept a few variants.
                lbl: Optional[str] = None
                w: Optional[float] = None
                if len(parts) == 4:
                    # Could be label or weight. Prefer weight if it parses cleanly.
                    try:
                        w = float(parts[3])
                    except ValueError:
                        lbl = parts[3]
                elif len(parts) >= 5:
                    lbl = parts[3]
                    try:
                        w = float(parts[4])
                    except ValueError:
                        w = None

                a, b = (u, v) if u <= v else (v, u)
                if lbl is not None:
                    edge_labels[(a, b)] = lbl
                if w is not None and float(w) != 1.0:
                    edge_weights[(a, b)] = float(w)
                continue

            # Unknown line type: ignore for robustness.
            continue

    # In case vertices were not explicitly listed, infer nodes from edges.
    if not nodes:
        s = set()
        for (u, v) in edges:
            s.add(u)
            s.add(v)
        nodes = sorted(s)

    return Graph(nodes=nodes, edges=edges, node_labels=node_labels or None, edge_labels=edge_labels or None, edge_weights=edge_weights or None)


def write_lg(graph: Graph, path: str | Path, directed: bool = False, include_weight: bool = False) -> None:
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
    edge_weights = graph.edge_weights or {}

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
            w = edge_weights.get((u_orig, v_orig), edge_weights.get((v_orig, u_orig), 1.0))

            if include_weight and float(w) != 1.0:
                # If there is no label, we still emit a placeholder label to keep parsing unambiguous.
                if lbl == "":
                    f.write(f"e {u} {v} _ {float(w)}\n")
                else:
                    f.write(f"e {u} {v} {lbl} {float(w)}\n")
            else:
                if lbl == "":
                    f.write(f"e {u} {v}\n")
                else:
                    f.write(f"e {u} {v} {lbl}\n")