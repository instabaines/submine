from __future__ import annotations

from typing import Any, Iterable, List, Optional

from ..core.graph import Graph
from ..core.result import SubgraphPattern

_DEFAULT_LABELS = ("__default_node__", "__default_edge__")


def _gspan_pattern_graph(rec: dict[str, Any], directed: bool, node_inv: Optional[dict] = None,
                         edge_inv: Optional[dict] = None) -> Graph:
    """Build a pattern Graph from one gSpan backend record, keeping its labels.

    The backend reports integer labels (as written to the gSpan file). ``node_inv`` /
    ``edge_inv`` map them back to the caller's original labels; without them (native
    gSpan input) the integer labels are kept. Undirected patterns list each edge in
    both directions; it is added once.
    """
    def back(inv, lbl):
        if inv is None:
            return lbl
        orig = inv.get(lbl, lbl)
        return None if orig in _DEFAULT_LABELS else orig

    nlabels = rec.get("node_labels") or []
    elabels = rec.get("edge_labels") or []
    g = Graph(directed=directed)
    for i, nid in enumerate(rec["nodes"]):
        g.add_node(nid, back(node_inv, nlabels[i]) if i < len(nlabels) else None)
    seen = set()
    for i, (u, v) in enumerate(rec["edges"]):
        key = (u, v) if directed else (min(u, v), max(u, v))
        if key in seen:
            continue
        seen.add(key)
        g.add_edge(u, v, label=back(edge_inv, elabels[i]) if i < len(elabels) else None)
    return g


def build_gspan_patterns(records: Iterable[dict[str, Any]], directed: bool, node_inv: Optional[dict] = None,
                         edge_inv: Optional[dict] = None) -> List[SubgraphPattern]:
    patterns: List[SubgraphPattern] = []
    for pid, rec in enumerate(records):
        pattern_graph = _gspan_pattern_graph(rec, directed, node_inv, edge_inv)
        patterns.append(
            SubgraphPattern(
                pid=pid,
                graph=pattern_graph,
                support=int(rec["support"]),
                frequency=None,
                occurrences=[],
                attributes={
                    "num_vertices": pattern_graph.number_of_nodes(),
                    "graph_ids": rec["graph_ids"],
                },
            )
        )
    return patterns


def build_sopagrami_patterns(records: Iterable[dict[str, Any]], directed: bool) -> List[SubgraphPattern]:
    patterns: List[SubgraphPattern] = []
    for pid, rec in enumerate(records):
        node_labels = list(rec["node_labels"])
        edges_raw = list(rec["edges"])
        support = int(rec["full_support"])
        key = rec["key"]

        # SoPaGraMi pattern node ids are 0..k-1
        k = len(node_labels)
        nodes = list(range(k))
        pat_edges: list[tuple[int, int]] = []
        edge_labels: dict[tuple[int, int], Any] = {}

        for (a, b, el, dirflag) in edges_raw:
            ai = int(a)
            bi = int(b)
            u, v = (ai, bi) if ai <= bi else (bi, ai)
            pat_edges.append((u, v))
            # the pattern Graph is undirected; direction is encoded in the label
            label = f"{el}->" if directed and dirflag == 1 else el
            edge_labels[(u, v)] = label

        node_label_map = {i: lbl for i, lbl in enumerate(node_labels)}
        pat_graph = Graph(
            nodes=nodes,
            edges=pat_edges,
            node_labels=node_label_map,
            edge_labels=edge_labels,
            directed=directed,
        )
        patterns.append(
            SubgraphPattern(
                pid=pid,
                graph=pat_graph,
                support=support,
                frequency=None,
                occurrences=[],
                attributes={"key": key, "k": k, "num_edges": len(pat_edges)},
            )
        )

    return patterns
