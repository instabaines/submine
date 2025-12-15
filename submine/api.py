from pathlib import Path
from typing import Iterable, Union, Sequence
import networkx as nx

from .core.graph import Graph
from . import get_mining_algorithm as get_algorithm
from .core.result import MiningResult, SubgraphPattern

GraphLike = Union[Graph, nx.Graph]
GraphSourceLike = Union[
    Graph,
    Iterable[Graph],
    Sequence[Graph],
    Path,
    str,
    # later: DB handles, etc.
]


def _normalize_graph_source(source: GraphSourceLike) -> Iterable[Graph]:
    # 1. Already an internal Graph → wrap in list
    if isinstance(source, Graph):
        return [source]

    # 2. Path / str → load from file
    if isinstance(source, (str, Path)):
        path = Path(source)
        # decide by extension or user config
        if path.suffix == ".gspan" or ".data" in str(source):
            from .io.gspan import read_gspan_dataset
            return read_gspan_dataset(path)
        elif path.suffix.lower() == ".lg":
            # Load a single lg graph (streaming reader)
            from .io.sopagrami import read_lg
            return [read_lg(path)]
        elif path.suffix in {".edgelist", ".txt"}:
            from .io.common import read_edgelist_dataset
            return read_edgelist_dataset(path)
        # etc.
        raise ValueError(f"Unsupported graph file format: {path.suffix}")

    # 3. Iterable of Graphs → pass through
    try:
        it = iter(source)  # type: ignore
    except TypeError:
        pass
    else:
        # could be list[Graph], generator, GraphSource, etc.
        # You may want to sanity check items, but can be lazy.
        return it

    raise TypeError(f"Cannot interpret {type(source)} as a graph source")


def mine_subgraphs(
    data: GraphSourceLike,
    algorithm: str,
    min_support: int,
    **algo_params,
) -> MiningResult:
    """High-level convenience function for users.

    `data` can be:
      - a single Graph
      - an iterable of Graphs
      - a path to a graph dataset on disk
    """
    AlgoCls = get_algorithm(algorithm)
    miner = AlgoCls(**algo_params)

    # Special-case: if the user supplies a native .lg file for SoPaGraMi,
    # pass it through directly (no re-parse, no re-write).
    if isinstance(data, (str, Path)):
        path = Path(data)
        if path.suffix.lower() == ".lg" and getattr(miner, "name", "").lower() == "sopagrami":
            if not hasattr(miner, "mine_lg"):
                raise RuntimeError("SoPaGraMi miner does not expose mine_lg(); update the wrapper.")
            return miner.mine_lg(path, min_support=min_support)

    graphs = _normalize_graph_source(data)
    return miner.mine(graphs, min_support=min_support)
