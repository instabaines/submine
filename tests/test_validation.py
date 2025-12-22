from __future__ import annotations

import pytest

from submine.api import mine_subgraphs
from submine.algorithms.gspan import GSpanMiner
from submine.core.graph import Graph
from submine.errors import ParameterValidationError


def test_api_min_support_validation_raises():
    g = Graph(nodes=[0, 1], edges=[(0, 1)])
    with pytest.raises(ParameterValidationError):
        mine_subgraphs([g], algorithm="gspan", min_support=0)


def test_gspan_miner_validates_vertex_bounds():
    with pytest.raises(ParameterValidationError):
        _ = GSpanMiner(min_vertices=0)

    with pytest.raises(ParameterValidationError):
        _ = GSpanMiner(min_vertices=3, max_vertices=2)
