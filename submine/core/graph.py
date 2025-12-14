"""Graph data structure used throughout the submine library.

The :class:`Graph` class represents a simple undirected graph where each
vertex and edge can carry an optional label. It is intentionally
lightweight and does not depend on external libraries such as
`networkx`, though conversion utilities are provided if networkx is
available. Internally the graph is stored as an adjacency dictionary
mapping node identifiers to a list of neighbouring nodes with edge
labels.

You can add nodes and edges explicitly using :meth:`add_node` and
:meth:`add_edge`. The graph is undirected by default – adding an edge
between ``u`` and ``v`` automatically registers the symmetric edge.

Example usage::

    >>> from submine.core.graph import Graph
    >>> g = Graph()
    >>> g.add_node(1, label='C')
    >>> g.add_node(2, label='O')
    >>> g.add_edge(1, 2, label='single')
    >>> list(g.iter_edges())
    [(1, 2, 'single')]

"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, Iterator, List, Optional, Tuple


from typing import Iterable, Protocol

class GraphSource(Protocol):
    """Anything that can yield Graph objects."""
    def __iter__(self) -> Iterable[Graph]:
        ...

@dataclass(eq=True, frozen=True)
class Node:
    """Internal representation of a graph node.

    Parameters
    ----------
    id: Any
        Unique identifier for the node.
    label: Optional[str]
        An optional label attached to the node. Labels are used by
        certain algorithms to distinguish node types.
    data: Dict[str, Any]
        Arbitrary metadata associated with this node.
    """

    id: Any
    label: Optional[str] = None
    data: Dict[str, Any] | None = None

    def __post_init__(self) -> None:
        # Ensure data is a dict even if None was provided.
        if self.data is None:
            object.__setattr__(self, "data", {})


class Graph:
    """A simple undirected labelled graph.

    Nodes are indexed by unique identifiers (any hashable type). Each
    node can optionally carry a label and arbitrary metadata. Edges
    connect unordered pairs of nodes and can also carry a label. Self
    loops are not supported.
    """

    def __init__(self, nodes, edges, node_labels=None, edge_labels=None):
        # nodes: iterable of node ids (can be any hashable)
        # edges: iterable of (u, v) or (u, v, label)
        self.nodes = list(nodes)
        self.edges = list(edges)
        self.node_labels = node_labels  # dict[node_id] -> label (int/str), or None
        self.edge_labels = edge_labels  # dict[(u, v)] -> label (int/str), or None

    def add_node(self, node_id: Any, label: Optional[str] = None, **data: Any) -> Node:
        """Add a new node to the graph.

        If a node with the given identifier already exists its label and
        data will be updated with the provided values. Extra keyword
        arguments are stored in the node's data dictionary.

        Parameters
        ----------
        node_id: Any
            Unique identifier for the node.
        label: Optional[str], default None
            A label for the node (e.g., atom type).
        **data: dict
            Arbitrary metadata to attach to the node.

        Returns
        -------
        Node
            The created or updated node.
        """
        if node_id in self._nodes:
            node = self._nodes[node_id]
            # Update existing node
            if label is not None:
                node = Node(node.id, label=label, data=node.data.copy())
            if data:
                node.data.update(data)
            self._nodes[node_id] = node
            return node
        # New node
        node = Node(node_id, label=label, data=data or None)
        self._nodes[node_id] = node
        self._adj.setdefault(node_id, [])
        return node

    def add_edge(self, u: Any, v: Any, label: Optional[str] = None) -> None:
        """Add an undirected edge between two nodes.

        Parameters
        ----------
        u, v: Any
            Identifiers of the nodes to connect. They must have been
            previously added via :meth:`add_node`.
        label: Optional[str], default None
            Label associated with this edge (e.g., bond type). The label
            will be stored for both directions.

        Raises
        ------
        KeyError
            If either ``u`` or ``v`` are not in the graph.
        ValueError
            If ``u`` equals ``v`` (self loops are not supported).
        """
        if u == v:
            raise ValueError("Self loops are not supported.")
        if u not in self._nodes or v not in self._nodes:
            raise KeyError(f"Both nodes must exist in the graph before adding an edge: {u}, {v}")
        # Ensure adjacency lists exist
        self._adj.setdefault(u, [])
        self._adj.setdefault(v, [])
        # Append in both directions
        self._adj[u].append((v, label))
        self._adj[v].append((u, label))

    def nodes(self) -> Iterable[Node]:
        """Return an iterable over all nodes in the graph."""
        return self._nodes.values()

    def iter_edges(self) -> Iterator[Tuple[Any, Any, Optional[str]]]:
        """Iterate over edges once.

        Yields
        ------
        Tuple[Any, Any, Optional[str]]
            A tuple (u, v, label) for each undirected edge in the graph.
            Each edge is yielded once with u < v according to the
            default Python ordering for the identifiers.
        """
        seen = set()
        for u, neighbours in self._adj.items():
            for v, label in neighbours:
                if (v, u) not in seen:
                    seen.add((u, v))
                    yield (u, v, label)

    def degree(self, node_id: Any) -> int:
        """Return the degree of a node."""
        return len(self._adj.get(node_id, []))

    def number_of_nodes(self) -> int:
        return len(self.nodes)

    def number_of_edges(self) -> int:
        return len(self.edges)

    def neighbours(self, node_id: Any) -> List[Tuple[Any, Optional[str]]]:
        """Return a list of neighbours and associated edge labels for a node."""
        return self._adj.get(node_id, [])

    def to_networkx(self) -> "networkx.Graph":
        """Convert this graph into a networkx.Graph.

        This convenience method requires `networkx` to be installed. If
        the import fails a RuntimeError is raised.

        Returns
        -------
        networkx.Graph
            A networkx representation of this graph.
        """
        try:
            import networkx as nx  # type: ignore
        except ImportError:
            raise RuntimeError("networkx is required for to_networkx()") from None
        G = nx.Graph()
        for node in self.nodes():
            G.add_node(node.id, label=node.label, **node.data)
        for u, v, label in self.iter_edges():
            G.add_edge(u, v, label=label)
        return G

    @classmethod
    def from_networkx(cls, G: "networkx.Graph") -> "Graph":
        """Create a :class:`Graph` from a networkx Graph.

        Node labels will be retrieved from the ``label`` attribute if
        present; other attributes are stored in the node's data.

        Parameters
        ----------
        G: networkx.Graph
            A networkx graph instance.

        Returns
        -------
        Graph
            A new Graph instance containing the same nodes and edges.
        """
        graph = cls()
        for n, attrs in G.nodes(data=True):
            label = attrs.pop("label", None)
            graph.add_node(n, label=label, **attrs)
        for u, v, attrs in G.edges(data=True):
            label = attrs.get("label")
            graph.add_edge(u, v, label=label)
        return graph

    def __repr__(self) -> str:
        return f"Graph(num_nodes={self.number_of_nodes()}, num_edges={self.number_of_edges()})"