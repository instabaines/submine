"""Built‑in datasets for submine.

This package provides convenience loaders for commonly used benchmark
datasets in frequent subgraph mining research. Datasets are returned
as lists of :class:`~submine.core.graph.Graph` objects. The available
datasets are documented in :mod:`submine.datasets.loaders`.
"""

from .loaders import get_dataset  # noqa: F401

__all__ = ["get_dataset"]