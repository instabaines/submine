"""Core graph container used throughout *submine*.

Design goals
------------
1) Keep a lightweight, dependency-free representation.
2) Preserve the existing public surface used by wrappers:
   - ``Graph.nodes`` : list of node ids
   - ``Graph.edges`` : list of (u, v)
   - ``Graph.node_labels`` : dict[node_id] -> label (optional)
   - ``Graph.edge_labels`` : dict[(u, v)] -> label (optional)
3) Add *optional* edge weights without breaking unweighted algorithms.

Weights are stored as ``Graph.edge_weights`` (dict[(u, v)] -> float). If a
weight is missing for an edge, it is treated as 1.0.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Hashable, Iterable, Iterator, List, Optional, Protocol, Tuple


class GraphSource(Protocol):
    """Anything that can yield Graph objects."""

    def __iter__(self) -> Iterable["Graph"]:  # pragma: no cover
        ...


@dataclass(eq=True, frozen=True)
class Node:
    id: Any
    label: Optional[Any] = None
    data: Dict[str, Any] | None = None

    def __post_init__(self) -> None:
        if self.data is None:
            object.__setattr__(self, "data", {})


class Graph:
    """Lightweight labeled (optionally weighted) undirected graph."""

    def __init__(
        self,
        nodes: Optional[Iterable[Hashable]] = None,
        edges: Optional[Iterable[Tuple[Hashable, Hashable]]] = None,
        node_labels: Optional[Dict[Hashable, Any]] = None,
        edge_labels: Optional[Dict[Tuple[Hashable, Hashable], Any]] = None,
        edge_weights: Optional[Dict[Tuple[Hashable, Hashable], float]] = None,
    ) -> None:
        self.nodes: List[Hashable] = list(nodes) if nodes is not None else []
        self.edges: List[Tuple[Hashable, Hashable]] = list(edges) if edges is not None else []
        self.node_labels: Optional[Dict[Hashable, Any]] = node_labels
        self.edge_labels: Optional[Dict[Tuple[Hashable, Hashable], Any]] = edge_labels
        self.edge_weights: Optional[Dict[Tuple[Hashable, Hashable], float]] = edge_weights

        # For incremental construction APIs.
        self._nodes: Dict[Hashable, Node] = {}
        self._adj: Dict[Hashable, List[Tuple[Hashable, Optional[Any], float]]] = {}

        if self.nodes or self.edges:
            # Seed internal indices for add_node/add_edge compatibility.
            for nid in self.nodes:
                lbl = self.node_labels.get(nid) if self.node_labels else None
                self._nodes[nid] = Node(nid, lbl, {})
                self._adj.setdefault(nid, [])
            for (u, v) in self.edges:
                lbl = None
                if self.edge_labels is not None:
                    lbl = self.edge_labels.get((u, v), self.edge_labels.get((v, u)))
                w = 1.0
                if self.edge_weights is not None:
                    w = float(self.edge_weights.get((u, v), self.edge_weights.get((v, u), 1.0)))
                self._adj.setdefault(u, []).append((v, lbl, w))
                self._adj.setdefault(v, []).append((u, lbl, w))

    @property
    def is_weighted(self) -> bool:
        if not self.edge_weights:
            return False
        return any(float(w) != 1.0 for w in self.edge_weights.values())

    def add_node(self, node_id: Hashable, label: Optional[Any] = None, **data: Any) -> Node:
        if node_id in self._nodes:
            n = self._nodes[node_id]
            lbl = n.label if label is None else label
            merged = dict(n.data)
            merged.update(data)
            n2 = Node(node_id, lbl, merged)
            self._nodes[node_id] = n2
            if node_id not in self.nodes:
                self.nodes.append(node_id)
            if self.node_labels is not None:
                self.node_labels[node_id] = lbl
            self._adj.setdefault(node_id, [])
            return n2

        n = Node(node_id, label, data or None)
        self._nodes[node_id] = n
        self._adj.setdefault(node_id, [])
        if node_id not in self.nodes:
            self.nodes.append(node_id)
        if label is not None:
            if self.node_labels is None:
                self.node_labels = {}
            self.node_labels[node_id] = label
        return n

    def add_edge(
        self,
        u: Hashable,
        v: Hashable,
        label: Optional[Any] = None,
        weight: float = 1.0,
    ) -> None:
        if u == v:
            raise ValueError("Self loops are not supported.")
        if u not in self._nodes:
            self.add_node(u)
        if v not in self._nodes:
            self.add_node(v)

        self.edges.append((u, v))
        self._adj.setdefault(u, []).append((v, label, float(weight)))
        self._adj.setdefault(v, []).append((u, label, float(weight)))

        if label is not None:
            if self.edge_labels is None:
                self.edge_labels = {}
            self.edge_labels[(u, v)] = label

        if float(weight) != 1.0:
            if self.edge_weights is None:
                self.edge_weights = {}
            self.edge_weights[(u, v)] = float(weight)

    def iter_edges(self) -> Iterator[Tuple[Hashable, Hashable, Optional[Any]]]:
        """Iterate edges once, yielding (u, v, label).

        This preserves the legacy shape expected by existing writers.
        """
        seen: set[Tuple[Hashable, Hashable]] = set()
        for (u, v) in self.edges:
            a, b = (u, v) if u <= v else (v, u)
            if (a, b) in seen:
                continue
            seen.add((a, b))
            lbl = None
            if self.edge_labels is not None:
                lbl = self.edge_labels.get((u, v), self.edge_labels.get((v, u)))
            yield (a, b, lbl)

    def iter_edges_with_weights(
        self,
    ) -> Iterator[Tuple[Hashable, Hashable, Optional[Any], float]]:
        """Iterate edges once, yielding (u, v, label, weight)."""
        seen: set[Tuple[Hashable, Hashable]] = set()
        for (u, v) in self.edges:
            a, b = (u, v) if u <= v else (v, u)
            if (a, b) in seen:
                continue
            seen.add((a, b))
            lbl = None
            if self.edge_labels is not None:
                lbl = self.edge_labels.get((u, v), self.edge_labels.get((v, u)))
            w = 1.0
            if self.edge_weights is not None:
                w = float(self.edge_weights.get((u, v), self.edge_weights.get((v, u), 1.0)))
            yield (a, b, lbl, w)

    def number_of_nodes(self) -> int:
        return len(self.nodes)

    def number_of_edges(self) -> int:
        return len(self.edges)

    def __repr__(self) -> str:
        return f"Graph(num_nodes={self.number_of_nodes()}, num_edges={self.number_of_edges()}, weighted={self.is_weighted})"
