from __future__ import annotations

from pathlib import Path
from typing import Iterable, Union, Sequence, Optional
import networkx as nx

import inspect
from .utils.checks import assert_regular_file
from .errors import SubmineInputError, ParameterValidationError


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


def _accepted_kwargs(callable_obj):
    sig = inspect.signature(callable_obj)
    accepted = set()
    has_var_kw = False
    for name, p in sig.parameters.items():
        if p.kind == inspect.Parameter.VAR_KEYWORD:
            has_var_kw = True
        elif p.kind in (inspect.Parameter.POSITIONAL_OR_KEYWORD, inspect.Parameter.KEYWORD_ONLY):
            accepted.add(name)
    return accepted, has_var_kw

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
    # Split kwargs between __init__ and mine()
    init_keys, init_var = _accepted_kwargs(AlgoCls.__init__)
    mine_keys, mine_var = _accepted_kwargs(AlgoCls.mine)

    init_params = {}
    run_params = {}
    unknown = {}

    for k, v in algo_params.items():
        in_init = (k in init_keys) or init_var
        in_mine = (k in mine_keys) or mine_var

        # Prefer explicit match if both accept (rare but possible)
        if k in init_keys and k in mine_keys:
            # Policy choice: treat as runtime override
            run_params[k] = v
        elif k in init_keys:
            init_params[k] = v
        elif k in mine_keys:
            run_params[k] = v
        else:
            unknown[k] = v

    if unknown:
        raise ParameterValidationError(
            f"Unsupported parameters for algorithm='{algorithm}': {sorted(unknown.keys())}"
        )
    miner = AlgoCls(**init_params)

    if not isinstance(min_support, int) or min_support <= 0:
        raise ParameterValidationError(f"min_support must be a positive integer; got {min_support!r}")

    # If user provided a path, and the miner declares a native on-disk format,
    # transcode directly to that format (only when needed) and call mine_native().
    if isinstance(data, (str, Path)):
        from .io.transcode import detect_format, transcode_path
        from .io.common import temporary_directory

        src_path = assert_regular_file(Path(data))
        src_fmt: Optional[str]
        try:
            src_fmt = detect_format(src_path)
        except Exception:
            src_fmt = None

        expected = getattr(miner, "expected_input_format", None)
        if expected is not None:
            if src_fmt == expected:
                return miner.mine_native(src_path, min_support=min_support, **run_params)

            # Not in the miner's native format: transcode once to native file.
            with temporary_directory() as tmp:
                suffix = ".lg" if expected == "lg" else ".data"
                native_path = tmp / f"native{suffix}"
                transcode_path(src_path, native_path, dst_fmt=expected, src_fmt=src_fmt)
                return miner.mine_native(native_path, min_support=min_support, **run_params)

    graphs = _normalize_graph_source(data)
    return miner.mine(graphs, min_support=min_support,**run_params)
