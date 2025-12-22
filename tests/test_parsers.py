from __future__ import annotations

from pathlib import Path

import pytest

from submine.core.graph import Graph
from submine.io.common import read_edgelist_dataset
from submine.io.gspan import read_gspan_dataset, write_gspan_dataset
from submine.io.sopagrami import read_lg
from submine.errors import ResourceLimitError


def data_path(name: str) -> Path:
    return Path(__file__).parent / "data" / name


def test_read_edgelist_dataset_triangle():
    graphs = read_edgelist_dataset(data_path("sample.edgelist"))
    assert len(graphs) == 1
    g = graphs[0]
    assert set(g.nodes) == {0, 1, 2}
    assert set(map(tuple, g.edges)) == {(0, 1), (1, 2), (2, 0)}


def test_read_gspan_dataset_basic():
    graphs = read_gspan_dataset(data_path("sample.data.x"))
    assert len(graphs) == 1
    g = graphs[0]
    assert set(g.nodes) == {0, 1}
    assert set(map(tuple, g.edges)) == {(0, 1)}
    assert g.node_labels is not None
    assert g.node_labels[0] == 2


def test_read_lg_parses_labels_and_weights():
    g = read_lg(data_path("sample.lg"))
    assert set(g.nodes) == {0, 1, 2}
    assert g.node_labels is not None
    assert g.node_labels[0] == "A"
    assert g.edge_labels is not None
    assert g.edge_labels[(0, 1)] == "rel"
    assert g.edge_weights is not None
    assert g.edge_weights[(0, 1)] == pytest.approx(0.5)


def test_write_gspan_dataset_is_deterministic(tmp_path: Path):
    g = Graph(
        nodes=[0, 1, 2],
        edges=[(0, 1), (1, 2)],
        node_labels={0: "X", 1: "X", 2: "Y"},
        edge_labels={(0, 1): "R", (1, 2): "R"},
    )
    out = tmp_path / "out.data"
    write_gspan_dataset([g], out)
    txt = out.read_text(encoding="utf-8")

    # Label remapping starts at 2; this exact output is stable for this input.
    assert txt.strip().splitlines() == [
        "t # 0",
        "v 0 2",
        "v 1 2",
        "v 2 3",
        "e 0 1 2",
        "e 1 2 2",
        "t # -1",
    ]


def test_streaming_limits_max_lines(monkeypatch, tmp_path: Path):
    # Force a very small max-lines cap.
    monkeypatch.setenv("SUBMINE_MAX_LINES", "2")

    p = tmp_path / "big.edgelist"
    p.write_text("0 1\n1 2\n2 3\n", encoding="utf-8")

    with pytest.raises(ResourceLimitError):
        _ = read_edgelist_dataset(p)
