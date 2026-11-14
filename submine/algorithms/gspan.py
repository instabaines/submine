from __future__ import annotations

from pathlib import Path
from typing import List, Optional

from ..utils.checks import safe_read_text, assert_regular_file
from ..errors import ParameterValidationError

import tempfile
import time

from .base import SubgraphMiner, register
from ..core.graph import Graph
from ..core.result import MiningResult, SubgraphPattern
from ..io.gspan import write_gspan_dataset, convert_gspan_graph, _build_label_maps
from typing import Iterable


_DEFAULT_LABELS = ("__default_node__", "__default_edge__")


def _pattern_graph(rec: dict, directed: bool, node_inv: Optional[dict] = None,
                   edge_inv: Optional[dict] = None) -> Graph:
    """Build a pattern Graph from one backend record, keeping its labels.

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
        lbl = back(node_inv, nlabels[i]) if i < len(nlabels) else None
        g.add_node(nid, lbl)
    seen = set()
    for i, (u, v) in enumerate(rec["edges"]):
        key = (u, v) if directed else (min(u, v), max(u, v))
        if key in seen:
            continue
        seen.add(key)
        lbl = back(edge_inv, elabels[i]) if i < len(elabels) else None
        g.add_edge(u, v, label=lbl)
    return g

@register
class GSpanMiner(SubgraphMiner):
    name = "gspan"
    expected_input_format = "gspan"

    def __init__(
        self,
        min_support: int = 2,
        directed: bool = False,
        min_vertices: int = 1,
        max_vertices: Optional[int] = None,
        visualize: bool = False,
        write_out: bool = True,
        verbose: bool = False,
    ) -> None:
        super().__init__(verbose=verbose)
        # Parameter validation (publish-safe defaults)
        if not isinstance(min_support, int) or min_support <= 0:
            raise ParameterValidationError(f"min_support must be a positive int; got {min_support!r}")
        if not isinstance(min_vertices, int) or min_vertices < 1:
            raise ParameterValidationError(f"min_vertices must be an int >= 1; got {min_vertices!r}")
        if max_vertices is not None:
            if not isinstance(max_vertices, int) or max_vertices < min_vertices:
                raise ParameterValidationError(
                    f"max_vertices must be None or an int >= min_vertices ({min_vertices}); got {max_vertices!r}"
                )
        self.min_support = min_support
        self.directed = directed
        self.min_vertices = min_vertices
        self.max_vertices = max_vertices
        self.visualize = visualize
        self.write_out = write_out


    def _run_on_dataset(self, db_path: Path, support: int):
        from . import gspan_cpp as gspan_mine

        t0 = time.time()

        db_path = assert_regular_file(db_path)
        data = safe_read_text(db_path)

        # TODO: plumb through additional kwargs once exposed by the binding.
        res = gspan_mine.mine_from_string(
            data,
            minsup=support,
            directed=self.directed,
            maxpat_min=self.min_vertices,
            maxpat_max=self.max_vertices if self.max_vertices is not None else 0xFFFFFFFF,
        )

        runtime = time.time() - t0
        return runtime, res

    def mine(self, graphs: List[Graph], min_support: Optional[int] = None, **kwargs) -> MiningResult:
        graphs = list(self._handle_weights(graphs))
        support = int(min_support if min_support is not None else self.min_support)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            db_path = tmpdir_path / "gspan_db.data"

            # write graphs in gspan format
            write_gspan_dataset(graphs, db_path)
            runtime, gs = self._run_on_dataset(db_path, support)

        # invert the label maps used by write_gspan_dataset, to report original labels
        node_map, edge_map = _build_label_maps(graphs)
        node_inv = {v: k for k, v in node_map.items()}
        edge_inv = {v: k for k, v in edge_map.items()}

        patterns = []
        for pid,rec in enumerate(gs):
            pattern_graph = _pattern_graph(rec, self.directed, node_inv, edge_inv)
            support = rec['support']
            patterns.append(
                SubgraphPattern(
                    pid=pid,
                    graph=pattern_graph,
                    support=support,
                    frequency=None,
                    occurrences=[],  # can fill later if track embeddings
                    attributes={
                        "num_vertices": pattern_graph.number_of_nodes(),
                        "graph_ids": rec["graph_ids"],
                    },
                )
            )
    

        return MiningResult(
            patterns=patterns,
            algorithm=self.name,
            params=dict(
                min_support=support,
                directed=self.directed,
                min_vertices=self.min_vertices,
                max_vertices=self.max_vertices,
                visualize=self.visualize,
                write_out=self.write_out,
            ),
            runtime=runtime,
            metadata={"backend": "gspan-mining"},
        )

    def mine_native(self, path: str | Path, min_support: int, **kwargs) -> MiningResult:
        """Run gSpan directly on a user-supplied gSpan dataset file."""
        db_path = Path(path)
        support = int(min_support)
        runtime, gs = self._run_on_dataset(db_path, support)

        patterns = []
        for pid,rec in enumerate(gs):
            pattern_graph = _pattern_graph(rec, self.directed)
            support = rec['support']
            patterns.append(
                SubgraphPattern(
                    pid=pid,
                    graph=pattern_graph,
                    support=support,
                    frequency=None,
                    occurrences=[],  # can fill later if track embeddings
                    attributes={
                        "num_vertices": pattern_graph.number_of_nodes(),
                        "graph_ids": rec["graph_ids"],
                    },
                )
            )

        return MiningResult(
            patterns=patterns,
            algorithm=self.name,
            params=dict(
                min_support=support,
                directed=self.directed,
                min_vertices=self.min_vertices,
                max_vertices=self.max_vertices,
                visualize=self.visualize,
                write_out=self.write_out,
                input_format="gspan",
            ),
            runtime=runtime,
            metadata={"backend": "gspan-mining", "input_dataset": str(db_path)},
        )