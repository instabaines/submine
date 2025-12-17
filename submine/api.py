from __future__ import annotations

from pathlib import Path
from typing import Iterable, Union, Sequence, Optional
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
        from .io.transcode import load_graphs

        return load_graphs(Path(source))

    # 3. Iterable of Graphs → pass through
    try:
        it = iter(source)  # type: ignore
    except TypeError:
        pass
    else:
        # could be list[Graph], generator, GraphSource, etc.
        #TODO:  sanity check items, but can be lazy.
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

    # If user provided a path, and the miner declares a native on-disk format,
    # transcode directly to that format (only when needed) and call mine_native().
    if isinstance(data, (str, Path)):
        from .io.transcode import detect_format, transcode_path
        from .io.common import temporary_directory

        src_path = Path(data)
        src_fmt: Optional[str]
        try:
            src_fmt = detect_format(src_path)
        except Exception:
            src_fmt = None

        expected = getattr(miner, "expected_input_format", None)
        if expected is not None:
            if src_fmt == expected:
                return miner.mine_native(src_path, min_support=min_support)

            # Not in the miner's native format: transcode once to native file.
            with temporary_directory() as tmp:
                suffix = ".lg" if expected == "lg" else ".data"
                native_path = tmp / f"native{suffix}"
                transcode_path(src_path, native_path, dst_fmt=expected, src_fmt=src_fmt)
                return miner.mine_native(native_path, min_support=min_support)

    graphs = _normalize_graph_source(data)
    return miner.mine(graphs, min_support=min_support)
