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
from typing import (
    Any,
    Dict,
    Hashable,
    Iterable,
    Iterator,
    List,
    Optional,
    Protocol,
    Tuple,
    Union,
)

class GraphSource(Protocol):
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

NodeSpec = Union[Hashable, Tuple[Hashable, Any], Node]
Edge2 = Tuple[Hashable, Hashable]
Edge3 = Tuple[Hashable, Hashable, Any]
Edge4 = Tuple[Hashable, Hashable, Any, float]
EdgeSpec = Union[Edge2, Edge3, Edge4]


class Graph:
    """Lightweight labeled (optionally weighted) graph.

    Ergonomic construction:
      - nodes can be: [id], [(id, label)], [Node(...)]
      - edges can be: [(u,v)], [(u,v,label)], [(u,v,label,weight)]

    Canonical storage:
      - For undirected graphs, labels/weights are stored under a canonical edge key
        so (u,v) and (v,u) are treated identically.
    """

    def __init__(
        self,
        nodes: Optional[Iterable[NodeSpec]] = None,
        edges: Optional[Iterable[EdgeSpec]] = None,
        node_labels: Optional[Dict[Hashable, Any]] = None,
        edge_labels: Optional[Dict[Tuple[Hashable, Hashable], Any]] = None,
        edge_weights: Optional[Dict[Tuple[Hashable, Hashable], float]] = None,
        directed: bool = False,
    ) -> None:
        self.directed = bool(directed)

        # Public-ish legacy fields (kept for compatibility).
        self.nodes: List[Hashable] = []
        self.edges: List[Tuple[Hashable, Hashable]] = []

        # Always keep these as dicts (None is what causes the “never updated” behavior).
        self.node_labels: Dict[Hashable, Any] = dict(node_labels) if node_labels else {}
        self.edge_labels: Dict[Tuple[Hashable, Hashable], Any] = dict(edge_labels) if edge_labels else {}
        self.edge_weights: Dict[Tuple[Hashable, Hashable], float] = dict(edge_weights) if edge_weights else {}

        # Internal indices for incremental APIs.
        self._nodes: Dict[Hashable, Node] = {}
        self._adj: Dict[Hashable, List[Tuple[Hashable, Optional[Any], float]]] = {}

        # 1) Ingest nodes (if provided)
        if nodes is not None:
            for spec in nodes:
                nid, lbl, data = self._parse_node_spec(spec)
                # Merge explicit node_labels mapping precedence:
                # - if node_labels contains nid, it wins
                if nid in self.node_labels:
                    lbl = self.node_labels[nid]
                self._seed_node(nid, lbl, data)

        # 2) Ingest edges (if provided), seeding missing nodes and labels/weights.
        if edges is not None:
            for espec in edges:
                u, v, elbl, w = self._parse_edge_spec(espec)

                # ensure endpoints exist
                if u not in self._nodes:
                    self._seed_node(u, self.node_labels.get(u), {})
                if v not in self._nodes:
                    self._seed_node(v, self.node_labels.get(v), {})

                self.edges.append((u, v))
                self._add_adj(u, v, elbl, w)

                # store label/weight canonically for undirected graphs
                if elbl is not None:
                    self.edge_labels[self._edge_key(u, v)] = elbl
                if float(w) != 1.0:
                    self.edge_weights[self._edge_key(u, v)] = float(w)

        # 3) Materialize legacy nodes list in stable order:
        #    - keep existing `self.nodes` insertion order as nodes were seeded.
        self.nodes = list(self._nodes.keys())

    # --------------------------
    # Parsing / canonicalization
    # --------------------------

    def _parse_node_spec(self, spec: NodeSpec) -> Tuple[Hashable, Optional[Any], Dict[str, Any]]:
        if isinstance(spec, Node):
            return spec.id, spec.label, dict(spec.data or {})
        if isinstance(spec, tuple) and len(spec) == 2:
            nid, lbl = spec
            return nid, lbl, {}
        return spec, None, {}

    def _parse_edge_spec(self, espec: EdgeSpec) -> Tuple[Hashable, Hashable, Optional[Any], float]:
        if len(espec) == 2:
            u, v = espec  # type: ignore[misc]
            return u, v, None, 1.0
        if len(espec) == 3:
            u, v, lbl = espec  # type: ignore[misc]
            return u, v, lbl, 1.0
        if len(espec) == 4:
            u, v, lbl, w = espec  # type: ignore[misc]
            return u, v, lbl, float(w)
        raise ValueError(f"Unsupported edge spec: {espec!r}")

    def _edge_key(self, u: Hashable, v: Hashable) -> Tuple[Hashable, Hashable]:
        """Canonical key for edge label/weight dictionaries."""
        if self.directed:
            return (u, v)
        # For undirected graphs, avoid relying on comparability (u <= v).
        # Use a stable ordering based on type name + repr as a fallback.
        su = (type(u).__name__, repr(u))
        sv = (type(v).__name__, repr(v))
        return (u, v) if su <= sv else (v, u)

    # --------------------------
    # Seeding helpers
    # --------------------------

    def _seed_node(self, node_id: Hashable, label: Optional[Any], data: Dict[str, Any]) -> None:
        if node_id in self._nodes:
            # merge: keep existing, update label/data if provided
            n = self._nodes[node_id]
            new_label = n.label if label is None else label
            merged = dict(n.data)
            merged.update(data or {})
            self._nodes[node_id] = Node(node_id, new_label, merged)
        else:
            self._nodes[node_id] = Node(node_id, label, data or None)

        self._adj.setdefault(node_id, [])

        if label is not None:
            self.node_labels[node_id] = label

    def _add_adj(self, u: Hashable, v: Hashable, label: Optional[Any], weight: float) -> None:
        self._adj.setdefault(u, []).append((v, label, float(weight)))
        if not self.directed:
            self._adj.setdefault(v, []).append((u, label, float(weight)))

    # --------------------------
    # Public API
    # --------------------------

    @property
    def is_weighted(self) -> bool:
        return any(float(w) != 1.0 for w in self.edge_weights.values())

    def add_node(self, node_id: Hashable, label: Optional[Any] = None, **data: Any) -> Node:
        self._seed_node(node_id, label, dict(data))
        if node_id not in self.nodes:
            self.nodes.append(node_id)
        return self._nodes[node_id]

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
        self._add_adj(u, v, label, float(weight))

        if label is not None:
            self.edge_labels[self._edge_key(u, v)] = label
        if float(weight) != 1.0:
            self.edge_weights[self._edge_key(u, v)] = float(weight)

    def iter_edges(self) -> Iterator[Tuple[Hashable, Hashable, Optional[Any]]]:
        """Iterate edges once, yielding (u, v, label)."""
        seen: set[Tuple[Hashable, Hashable]] = set()
        for (u, v) in self.edges:
            k = self._edge_key(u, v)
            if not self.directed:
                if k in seen:
                    continue
                seen.add(k)
                a, b = k
            else:
                a, b = u, v

            lbl = self.edge_labels.get(self._edge_key(u, v))
            yield (a, b, lbl)

    def iter_edges_with_weights(self) -> Iterator[Tuple[Hashable, Hashable, Optional[Any], float]]:
        """Iterate edges once, yielding (u, v, label, weight)."""
        seen: set[Tuple[Hashable, Hashable]] = set()
        for (u, v) in self.edges:
            k = self._edge_key(u, v)
            if not self.directed:
                if k in seen:
                    continue
                seen.add(k)
                a, b = k
            else:
                a, b = u, v

            key = self._edge_key(u, v)
            lbl = self.edge_labels.get(key)
            w = float(self.edge_weights.get(key, 1.0))
            yield (a, b, lbl, w)

    def number_of_nodes(self) -> int:
        return len(self.nodes)

    def number_of_edges(self) -> int:
        return len(self.edges)

    def __repr__(self) -> str:
        return (
            f"Graph(num_nodes={self.number_of_nodes()}, "
            f"num_edges={self.number_of_edges()}, "
            f"directed={self.directed}, weighted={self.is_weighted})"
        )
