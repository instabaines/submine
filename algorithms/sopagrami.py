"""Placeholder for SoPaGraMi algorithm wrapper.

SoPaGraMi (SOcial PAtterns Graph Miner) is an algorithm designed to
mine frequent subgraphs in dynamic or social network settings. Like
GraMi, the reference implementation is typically provided in another
language such as Java. This module defines a stub wrapper that
illustrates how such an algorithm might be integrated into the
submine framework. No actual mining is performed.
"""

from __future__ import annotations

from typing import List

from .base import SubgraphMiner, register
from ..core.graph import Graph
from ..core.result import MiningResult


@register
class SoPaGraMiMiner(SubgraphMiner):
    """Stub wrapper for the SoPaGraMi algorithm."""

    name = "sopagrami"

    def __init__(self, verbose: bool = False) -> None:
        super().__init__(verbose=verbose)

    def check_availability(self) -> None:
        # In a real implementation this would check for the presence of the
        # SoPaGraMi binary or jar.
        return None

    def mine(self, graphs: List[Graph], min_support: int, **kwargs) -> MiningResult:
        raise NotImplementedError(
            "SoPaGraMi integration is not yet implemented. Please add implementation."
        )