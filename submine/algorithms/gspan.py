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
from ..io.gspan import write_gspan_dataset,convert_gspan_graph
from typing import Iterable

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

        patterns = []
        for pid,rec in enumerate(gs):
            pattern_graph = Graph(edges=rec["edges"],nodes=rec['nodes'])
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
            pattern_graph = Graph(edges=rec["edges"],nodes=rec['nodes'])
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