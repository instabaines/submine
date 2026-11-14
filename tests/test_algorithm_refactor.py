from __future__ import annotations

from submine.algorithms.gspan import GSpanMiner
from submine.algorithms.pattern_utils import build_sopagrami_patterns
from submine.core.graph import Graph


def test_gspan_reports_requested_min_support(monkeypatch):
    miner = GSpanMiner(min_support=2)
    graphs = [Graph(nodes=[0, 1], edges=[(0, 1)])]

    def fake_run(_db_path, support):
        assert support == 3
        return 0.01, [
            {"nodes": [0, 1], "edges": [(0, 1)], "support": 9, "graph_ids": [0]},
            {"nodes": [0, 1, 2], "edges": [(0, 1), (1, 2)], "support": 7, "graph_ids": [0]},
        ]

    monkeypatch.setattr(miner, "_run_on_dataset", fake_run)

    result = miner.mine(graphs, min_support=3)

    assert result.params["min_support"] == 3
    assert [p.support for p in result.patterns] == [9, 7]


def test_sopagrami_pattern_builder_handles_directed_edge_marker():
    records = [
        {
            "node_labels": ["A", "B"],
            "edges": [(0, 1, "x", 1)],
            "full_support": 4,
            "key": "k1",
        }
    ]

    patterns = build_sopagrami_patterns(records, directed=True)
    graph = patterns[0].graph

    assert patterns[0].support == 4
    assert graph.edge_labels[(0, 1)] == "x->"


def test_weight_ignore_keeps_graph_consistent():
    miner = GSpanMiner()
    g = Graph(nodes=[0, 1], edges=[(0, 1)], edge_weights={(0, 1): 2.5})

    [unweighted] = list(miner._handle_weights([g]))

    assert unweighted.edge_weights == {}
    assert unweighted.is_weighted is False
