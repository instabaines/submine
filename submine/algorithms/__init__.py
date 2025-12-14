"""Algorithm implementations for submine.

Each submodule in this package implements a specific subgraph mining
algorithm. Modules are expected to define a subclass of
:class:`~submine.algorithms.base.SubgraphMiner` and register it via
:func:`~submine.algorithms.base.register`. Registered algorithms will
automatically appear in :func:`submine.get_algorithm`.

To avoid the cost of importing heavy dependencies at module import
time, algorithm modules should not perform expensive setup at the top
level. Instead they should defer initialization to the constructor or
:meth:`SubgraphMiner.check_availability`.

"""

from .base import SubgraphMiner  # noqa: F401

# Import algorithm modules so they can register themselves when this
# package is imported. Additional algorithms can be added here.
from .gspan import GSpanMiner  # noqa: F401
from .sopagrami import SoPaGraMiMiner  # noqa: F401

__all__ = ["SubgraphMiner", "GSpanMiner", "SoPaGraMiMiner"]