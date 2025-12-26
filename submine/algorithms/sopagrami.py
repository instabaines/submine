# submine/algorithms/sopagrami.py
from __future__ import annotations

import tempfile
import time
from pathlib import Path
from typing import Iterable, List, Optional

from .base import SubgraphMiner, register
from ..core.graph import Graph
from ..core.result import MiningResult, SubgraphPattern
from ..io.sopagrami import read_lg, write_lg
from ..errors import ParameterValidationError




@register
class SoPaGraMiMiner(SubgraphMiner):
    """
    Python wrapper around the C++ SoPaGraMi implementation.

    Note: SoPaGraMi mines frequent subgraphs from a *single* large graph,
    not a dataset of many graphs.
    """
    name = "sopagrami"
    expected_input_format = "lg"
    multi_graph_policy = "reject"

    def __init__(
        self,
        tau: int = 2,
        directed: bool = False,
        sorted_seeds: bool = True,
        num_threads: int = 0,
        compute_full_support: bool = True,
        verbose: bool = False,
    ) -> None:
        super().__init__(verbose=verbose)
        # Parameter validation (publish-safe defaults)
        if not isinstance(tau, int) or tau < 1:
            raise ParameterValidationError(f"tau must be an int >= 1; got {tau!r}")
        if not isinstance(num_threads, int) or num_threads < 0:
            raise ParameterValidationError(f"num_threads must be an int >= 0; got {num_threads!r}")
        self.tau = tau
        self.directed = directed
        self.sorted_seeds = sorted_seeds
        self.num_threads = num_threads
        self.compute_full_support = compute_full_support

    def check_availability(self):
        try:
            from . import sopagrami_cpp
        except ImportError as e:
            raise RuntimeError("SoPaGraMi backend not available") from e

    def mine(
        self,
        graphs: Iterable[Graph],
        min_support: Optional[int] = None,
        out_dir:str=None,dump_images_csv:bool=False,
                           max_images_per_vertex:int=50,dump_sample_embeddings:bool=False
    ) -> MiningResult:
        self.check_availability()

        # Handle weights explicitly (SoPaGraMi backend treats graphs as labeled, not weighted).
        graphs = self._handle_weights(graphs)
        # SoPaGraMi expects a single graph
        graphs_list = list(graphs)
        if len(graphs_list) != 1:
            raise ValueError(
                "SoPaGraMiMiner currently expects exactly one Graph (single large graph). "
                f"Got {len(graphs_list)}."
            )
        G = graphs_list[0]

        tau = int(min_support if min_support is not None else self.tau)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            lg_path = tmpdir_path / "graph.lg"

            # 1) write graph as .lg
            write_lg(G, lg_path, directed=self.directed)

            # 2) call C++ binding
            runtime, patterns_raw = self._run_backend_on_lg(lg_path, tau=tau,out_dir=out_dir,dump_images_csv=dump_images_csv,
                           max_images_per_vertex=max_images_per_vertex,dump_sample_embeddings=dump_sample_embeddings)

        # 3) Convert to our SubgraphPattern representation
        patterns: List[SubgraphPattern] = []
        for pid, pd in enumerate(patterns_raw):
            node_labels = list(pd["node_labels"])
            edges_raw = list(pd["edges"])
            support = int(pd["full_support"])
            key = pd["key"]

            # SoPaGraMi pattern node IDs are 0..k-1
            k = len(node_labels)
            nodes = list(range(k))

            # Build our Graph for the pattern
            pat_edges = []
            edge_labels = {}
            for (a, b, el, dirflag) in edges_raw:
                a = int(a)
                b = int(b)
                # Our Graph is undirected; we store the undirected edge,
                # and put direction info into the label if needed.
                u, v = (a, b) if a <= b else (b, a)
                pat_edges.append((u, v))

                label = el
                if self.directed and dirflag == 1:
                    # encode direction in the label for now
                    label = f"{el}->"
                edge_labels[(u, v)] = label

            node_label_map = {i: lbl for i, lbl in enumerate(node_labels)}

            pat_graph = Graph(
                nodes=nodes,
                edges=pat_edges,
                node_labels=node_label_map,
                edge_labels=edge_labels,
            )

            patterns.append(
                SubgraphPattern(
                    pid=pid,
                    graph=pat_graph,
                    support=support,
                    frequency=None,
                    occurrences=[],
                    attributes={
                        "key": key,
                        "k": k,
                        "num_edges": len(pat_edges),
                    },
                )
            )

        return MiningResult(
            patterns=patterns,
            algorithm=self.name,
            params={
                "tau": tau,
                "directed": self.directed,
                "sorted_seeds": self.sorted_seeds,
                "num_threads": self.num_threads,
                "compute_full_support": self.compute_full_support,
                

            },
            runtime=runtime,
            metadata={"backend": "sopagrami_cpp"},
        )

    def mine_native(self, lg_path: str | Path, min_support: Optional[int] = None, out_dir:str=None,dump_images_csv:bool=False,
                           max_images_per_vertex:int=50,dump_sample_embeddings:bool=False) -> MiningResult:
        """Run SoPaGraMi directly on a user-supplied ``.lg`` file.

        This avoids re-parsing/re-writing the file, which is important for
        large graphs and for preserving any optional attributes present in the
        original ``.lg``.
        """
        self.check_availability()
        lg_path = Path(lg_path)
        if lg_path.suffix.lower() != ".lg":
            raise ValueError(f"Expected a .lg file for SoPaGraMi; got: {lg_path}")

        tau = int(min_support if min_support is not None else self.tau)
        runtime, patterns_raw = self._run_backend_on_lg(lg_path, tau=tau,out_dir=out_dir,dump_images_csv=dump_images_csv,
                           max_images_per_vertex=max_images_per_vertex,dump_sample_embeddings=dump_sample_embeddings)

        # Convert patterns (same as in mine())
        patterns: List[SubgraphPattern] = []
        for pid, pd in enumerate(patterns_raw):
            node_labels = list(pd["node_labels"])
            edges_raw = list(pd["edges"])
            support = int(pd["full_support"])
            key = pd["key"]

            k = len(node_labels)
            nodes = list(range(k))

            pat_edges = []
            edge_labels = {}
            for (a, b, el, dirflag) in edges_raw:
                a = int(a)
                b = int(b)
                u, v = (a, b) if a <= b else (b, a)
                pat_edges.append((u, v))
                label = el
                if self.directed and dirflag == 1:
                    label = f"{el}->"
                edge_labels[(u, v)] = label

            node_label_map = {i: lbl for i, lbl in enumerate(node_labels)}
            pat_graph = Graph(nodes=nodes, edges=pat_edges, node_labels=node_label_map, edge_labels=edge_labels)

            patterns.append(
                SubgraphPattern(
                    pid=pid,
                    graph=pat_graph,
                    support=support,
                    frequency=None,
                    occurrences=[],
                    attributes={"key": key, "k": k, "num_edges": len(pat_edges)},
                )
            )

        return MiningResult(
            patterns=patterns,
            algorithm=self.name,
            params={
                "tau": tau,
                "directed": self.directed,
                "sorted_seeds": self.sorted_seeds,
                "num_threads": self.num_threads,
                "compute_full_support": self.compute_full_support,
                "input_format": "lg"
            },
            runtime=runtime,
            metadata={"backend": "sopagrami_cpp", "input_lg": str(lg_path)},
        )

    def _run_backend_on_lg(self, lg_path: Path, tau: int,out_dir:str=None,dump_images_csv:bool=False,
                           max_images_per_vertex:int=50,dump_sample_embeddings:bool=False):
        from . import sopagrami_cpp
        t0 = time.time()
        self.logger.debug("Running SoPaGraMi on %s", lg_path)
        if out_dir is None:
            out_dir = "sopagrami_result"
        patterns_raw = sopagrami_cpp.run_on_lg_file(
            str(lg_path),
            tau=tau,
            directed=self.directed,
            sorted_seeds=self.sorted_seeds,
            num_threads=self.num_threads,
            compute_full_support=self.compute_full_support,
            dump_images_csv = dump_images_csv,
            out_dir = out_dir,
            max_images_per_vertex = max_images_per_vertex,
            dump_sample_embeddings=dump_sample_embeddings



        )
        return time.time() - t0, patterns_raw