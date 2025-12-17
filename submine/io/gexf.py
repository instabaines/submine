"""GEXF reader.

We use NetworkX's GEXF parser to support common graph exports.

Mapping rules
-------------
- Node labels: if a node attribute named ``label`` exists, it is used; otherwise
  we use the node id as its label.
- Edge labels: if an edge attribute named ``label`` exists, it is used.
- Edge weights: if an edge attribute named ``weight`` exists and is not 1.0, it
  is stored in :attr:`submine.core.graph.Graph.edge_weights`.

Notes
-----
GEXF can represent directed graphs, multi-edges, and parallel edges. The internal
``Graph`` container in *submine* is currently undirected and does not preserve
parallel edges. We therefore:
  - coerce to an undirected simple graph
  - keep the first-seen label/weight for each undirected edge
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Hashable, Tuple

import networkx as nx

from ..core.graph import Graph


def read_gexf(path: str | Path) -> Graph:
    p = Path(path)
    g_nx = nx.read_gexf(p)

    # Coerce to undirected simple graph to match our internal container.
    if isinstance(g_nx, (nx.MultiGraph, nx.MultiDiGraph)):
        g_simple = nx.Graph()
        for u, v, data in g_nx.edges(data=True):
            if u == v:
                continue
            if g_simple.has_edge(u, v):
                continue
            g_simple.add_edge(u, v, **(data or {}))
        for n, data in g_nx.nodes(data=True):
            g_simple.add_node(n, **(data or {}))
        g_nx = g_simple
    else:
        g_nx = nx.Graph(g_nx)

    nodes = list(g_nx.nodes())

    node_labels: Dict[Hashable, Any] = {}
    for n, data in g_nx.nodes(data=True):
        if data is None:
            node_labels[n] = str(n)
        else:
            node_labels[n] = data.get("label", str(n))

    edges = []
    edge_labels: Dict[Tuple[Hashable, Hashable], Any] = {}
    edge_weights: Dict[Tuple[Hashable, Hashable], float] = {}

    for u, v, data in g_nx.edges(data=True):
        if u == v:
            continue
        a, b = (u, v) if str(u) <= str(v) else (v, u)
        key = (a, b)
        edges.append(key)

        if data:
            if "label" in data and data["label"] is not None:
                edge_labels[key] = data["label"]
            if "weight" in data and data["weight"] is not None:
                try:
                    w = float(data["weight"])
                except Exception:
                    w = 1.0
                if w != 1.0:
                    edge_weights[key] = w

    return Graph(
        nodes=nodes,
        edges=edges,
        node_labels=node_labels or None,
        edge_labels=edge_labels or None,
        edge_weights=edge_weights or None,
    )
