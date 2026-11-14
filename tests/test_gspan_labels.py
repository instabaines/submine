"""Edge/node labels survive the gSpan round trip (writer -> backend -> MiningResult)."""
from __future__ import annotations

from pathlib import Path

import pytest

from submine.core.graph import Graph
from submine.io.gspan import _build_label_maps, write_gspan_dataset


def _edge_lines(path: Path) -> list[str]:
    return [l for l in path.read_text().splitlines() if l.startswith("e ")]


def test_writer_uses_canonical_edge_key(tmp_path: Path):
    # Undirected keys are ordered by repr, so (5, 24) is stored under (24, 5).
    g = Graph(nodes=[(5, "a"), (24, "b")], edges=[(5, 24, "x")])
    out = tmp_path / "db.data"
    write_gspan_dataset([g], out)
    _, edge_map = _build_label_maps([g])
    assert _edge_lines(out) == [f"e 0 1 {edge_map['x']}"]


def test_writer_handles_partially_labeled_graphs(tmp_path: Path):
    g = Graph(nodes=[(0, "a"), 1, (2, "b")], edges=[(0, 1, "x"), (1, 2)])
    out = tmp_path / "db.data"
    write_gspan_dataset([g], out)
    node_map, edge_map = _build_label_maps([g])
    assert "__default_node__" in node_map and "__default_edge__" in edge_map
    assert len(_edge_lines(out)) == 2
    assert len(set(node_map.values())) == len(node_map)   # default does not collide with a real label
    assert len(set(edge_map.values())) == len(edge_map)


def test_mined_patterns_keep_labels():
    pytest.importorskip("submine.algorithms.gspan_cpp")
    from submine.api import mine_subgraphs

    graphs = [Graph(nodes=[(0, "A"), (1, "B"), (2, "C")], edges=[(0, 1, "x"), (1, 2, "y")]) for _ in range(3)]
    res = mine_subgraphs(graphs, algorithm="gspan", min_support=3)
    three = [p for p in res.patterns if p.graph.number_of_nodes() == 3]
    assert three, "expected the full 3-node pattern"
    g = three[0].graph
    assert sorted(g.node_labels.values()) == ["A", "B", "C"]
    assert sorted(g.edge_labels.values()) == ["x", "y"]
    assert len(g.edges) == 2          # undirected edges reported once
