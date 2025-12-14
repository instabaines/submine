"""Core data structures for the submine library.

This subpackage contains fundamental classes used throughout the
library. Currently the primary exported objects are
:class:`~submine.core.graph.Graph` for representing graphs and
:class:`~submine.core.result.FrequentSubgraph` for storing mining results.
"""

from .graph import Graph 
from .result import SubgraphPattern, MiningResult  

__all__ = ["Graph", "SubgraphPattern", "MiningResult"]