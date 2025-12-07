# submine/core/result.py
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from .graph import Graph

__all__ = [
    "SubgraphOccurrence",
    "SubgraphPattern",
    "MiningResult",
]


@dataclass
class SubgraphOccurrence:
    """
    One occurrence (embedding) of a pattern inside a particular input graph.

    Attributes
    ----------
    graph_id : int
        Index of the input graph in the dataset (0-based).
    node_mapping : Mapping[int, int]
        Mapping from pattern-local node ids -> node ids in the host graph.
        Pattern-local ids should be 0..k-1 for a k-node pattern.
    extra : dict
        Optional algorithm-specific information
        (e.g. edge mapping, score, position, etc.).
    """
    graph_id: int
    node_mapping: Mapping[int, int]
    extra: Dict[str, Any] = field(default_factory=dict)


@dataclass
class SubgraphPattern:
    """
    A mined subgraph pattern.

    Attributes
    ----------
    pid : int
        Pattern identifier (unique within a MiningResult).
    graph : Graph
        The pattern itself as a Graph object (pattern-local node ids 0..k-1).
    support : int
        Support count (e.g. number of graphs / embeddings where it appears).
    frequency : Optional[float]
        Relative frequency (e.g. support / num_graphs); optional.
    occurrences : list[SubgraphOccurrence]
        Concrete occurrences of this pattern in the input graphs (optional;
        some algorithms may not provide embeddings).
    attributes : dict
        Algorithm-specific attributes (e.g. score, interestingness, DFScode).
    """
    pid: int
    graph: Graph
    support: int
    frequency: Optional[float] = None
    occurrences: List[SubgraphOccurrence] = field(default_factory=list)
    attributes: Dict[str, Any] = field(default_factory=dict)

    def add_occurrence(self, occ: SubgraphOccurrence) -> None:
        self.occurrences.append(occ)


@dataclass
class MiningResult:
    """
    Container for the output of a subgraph mining run.

    Attributes
    ----------
    patterns : list[SubgraphPattern]
        List of mined patterns.
    algorithm : str
        Name of the algorithm that produced these patterns.
    params : dict
        Parameters used for the run (min_support, max_size, etc.).
    runtime : Optional[float]
        Wall-clock runtime in seconds, if measured.
    metadata : dict
        Additional metadata (dataset info, version, logs, etc.).
    """
    patterns: List[SubgraphPattern]
    algorithm: str
    params: Dict[str, Any] = field(default_factory=dict)
    runtime: Optional[float] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.patterns)

    def top_k(self, k: int, key: str = "support") -> List[SubgraphPattern]:
        """
        Return the top-k patterns sorted by a given key.

        Parameters
        ----------
        k : int
            Number of patterns to return.
        key : str
            Sort key: 'support', 'frequency', or any numeric attribute name
            stored in pattern.attributes[key].

        Returns
        -------
        list[SubgraphPattern]
        """
        def key_fn(p: SubgraphPattern):
            if key == "support":
                return p.support
            if key == "frequency":
                return p.frequency if p.frequency is not None else 0.0
            # fallback to attributes
            val = p.attributes.get(key, 0.0)
            return float(val)

        return sorted(self.patterns, key=key_fn, reverse=True)[:k]
