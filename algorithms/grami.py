"""Wrapper for the GraMi algorithm.

GraMi (Graph Miner) is a Java‑based algorithm for frequent subgraph
mining. This module defines a placeholder wrapper that can be extended
to call the GraMi implementation externally. To use GraMi you need to
download the ``grami.jar`` and ensure that a Java runtime is
available. Once the jar is placed in your environment you can extend
this class to invoke it via :func:`~submine.algorithms.base.SubgraphMiner.run_external`.

The current implementation is a stub and simply raises
:class:`NotImplementedError` when mining is attempted.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

from .base import SubgraphMiner, register
from ..core.graph import Graph
from ..core.result import MiningResult


@register
class GraMiMiner(SubgraphMiner):
    """Stub wrapper for the GraMi frequent subgraph miner."""

    name = "grami"

    def __init__(self, jar_path: str | None = None, verbose: bool = False) -> None:
        super().__init__(verbose=verbose)
        # Path to the GraMi jar. If not provided it will be looked up
        # based on an environment variable or default location.
        self.jar_path = Path(jar_path) if jar_path else None

    def check_availability(self) -> None:
        # Example check: verify that Java is installed and the jar is present
        self.logger.debug("Checking availability of GraMi and Java runtime")
        # For demonstration we do not implement actual checks here
        return None

    def mine(self, graphs: List[Graph], min_support: int, **kwargs) -> MiningResult:
        """Run the GraMi algorithm on the provided graphs.

        Raises
        ------
        NotImplementedError
            Always, since this is a stub implementation.
        """
        raise NotImplementedError(
            "GraMi integration is not yet implemented. Provide a jar and implement the invocation."
        )