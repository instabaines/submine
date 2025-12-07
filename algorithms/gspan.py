from __future__ import annotations

from pathlib import Path
from typing import List, Optional
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
        self.min_support = min_support
        self.directed = directed
        self.min_vertices = min_vertices
        self.max_vertices = max_vertices
        self.visualize = visualize
        self.write_out = write_out
    
        try:
            from  gSpan.gspan_mining.config import parser  # type: ignore
            from  gSpan.gspan_mining.main import main as gspan_main  # type: ignore
        except ImportError as e:
            raise ImportError(
                "GSpanMiner requires the 'gspan-ming' package.\n"
                "Install it with: pip install gspan-mining"
            ) from e

        self._parser = parser
        self._gspan_main = gspan_main

    def mine(self, graphs: List[Graph], min_support: Optional[int] = None, **kwargs) -> MiningResult:
        support = int(min_support if min_support is not None else self.min_support)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            db_path = tmpdir_path / "gspan_db"

            # write graphs in gspan format
            write_gspan_dataset(graphs, db_path)

            args = [
                "-s", str(support),
                "-d", str(self.directed),
                "-l", str(self.min_vertices),
                "-p", str(self.visualize),
                "-w", str(self.write_out),
                str(db_path),
            ]
            if self.max_vertices is not None:
                args.extend(["-u", str(self.max_vertices)])

            FLAGS, _ = self._parser.parse_known_args(args=args)
            t0 = time.time()
            self.logger.debug("Calling gSpan with args: %s", " ".join(args))
            gs = self._gspan_main(FLAGS)  # noqa: F841  # we'll hook patterns later
            runtime = time.time() - t0

        
        patterns = []
        for pid, rec in enumerate(gs.results):
            g_span_graph = rec["graph"]          # gspan.graph.Graph
            support = rec["support"]
            pattern_graph = convert_gspan_graph(g_span_graph)

            patterns.append(
                SubgraphPattern(
                    pid=pid,
                    graph=pattern_graph,
                    support=support,
                    frequency=None,
                    occurrences=[],  # can fill later if you track embeddings
                    attributes={
                        "num_vertices": rec["num_vertices"],
                        "description": rec["description"],
                        "where": rec["where"],
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
