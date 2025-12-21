"""Dataset loaders for submine.

This module contains functions to load well‑known benchmark datasets
such as MUTAG, ENZYMES and more. Where possible it attempts to use
existing libraries to fetch the datasets (e.g. torch_geometric or
networkx), but falls back to a simple synthetic dataset if these are
unavailable. To add a new dataset, implement a function named
``load_<datasetname>()`` that returns a list of
:class:`~submine.core.graph.Graph` objects and register it in the
``_DATASETS`` dictionary.
"""

from __future__ import annotations

from typing import Callable, Dict, List

from ..core.graph import Graph


def _load_toy() -> List[Graph]:
    """Return a small synthetic dataset useful for testing.

    The dataset contains two graphs: a triangle and a path of length 3.
    Node labels are single characters.
    """
    # Graph 1: triangle A-B-C-A
    g1 = Graph()
    g1.add_node(0, label="A")
    g1.add_node(1, label="B")
    g1.add_node(2, label="C")
    g1.add_edge(0, 1, label="ab")
    g1.add_edge(1, 2, label="bc")
    g1.add_edge(2, 0, label="ca")
    # Graph 2: path D-E-F
    g2 = Graph()
    g2.add_node(0, label="D")
    g2.add_node(1, label="E")
    g2.add_node(2, label="F")
    g2.add_edge(0, 1, label="de")
    g2.add_edge(1, 2, label="ef")
    return [g1, g2]


def _load_mutag() -> List[Graph]:
    """Load the MUTAG dataset if torch_geometric is available.

    If the dataset cannot be downloaded or the library is missing a
    NotImplementedError is raised.
    """
    try:
        from torch_geometric.datasets import TUDataset  # type: ignore
        import torch_geometric.utils  # noqa: F401  # type: ignore
    except Exception as e:
        raise NotImplementedError("Loading MUTAG requires torch_geometric and internet access") from e
    # Fetch dataset
    dataset = TUDataset(root="/tmp/MUTAG", name="MUTAG")
    graphs: List[Graph] = []
    for data in dataset:
        # Convert to networkx graph then to our Graph
        import networkx as nx  # type: ignore
        g_nx = nx.Graph()
        # Nodes with labels from x attribute (one-hot) if present
        for i in range(data.num_nodes):
            label = None
            if hasattr(data, "node_label") and data.node_label is not None:
                label = str(int(data.node_label[i]))
            g_nx.add_node(i, label=label)
        # Edges
        edge_attr = data.edge_attr if data.edge_attr is not None else None
        for j in range(data.edge_index.size(1)):
            u = int(data.edge_index[0, j])
            v = int(data.edge_index[1, j])
            label = None
            if edge_attr is not None:
                label = str(int(edge_attr[j].item()))
            g_nx.add_edge(u, v, label=label)
        graphs.append(Graph.from_networkx(g_nx))
    return graphs


def _load_enzymes() -> List[Graph]:
    """Load the ENZYMES dataset if torch_geometric is available."""
    try:
        from torch_geometric.datasets import TUDataset  # type: ignore
        import torch_geometric.utils  # noqa: F401  # type: ignore
    except Exception as e:
        raise NotImplementedError("Loading ENZYMES requires torch_geometric and internet access") from e
    dataset = TUDataset(root="/tmp/ENZYMES", name="ENZYMES")
    graphs: List[Graph] = []
    for data in dataset:
        import networkx as nx  # type: ignore
        g_nx = nx.Graph()
        for i in range(data.num_nodes):
            label = None
            if hasattr(data, "node_label") and data.node_label is not None:
                label = str(int(data.node_label[i]))
            g_nx.add_node(i, label=label)
        edge_attr = data.edge_attr if data.edge_attr is not None else None
        for j in range(data.edge_index.size(1)):
            u = int(data.edge_index[0, j])
            v = int(data.edge_index[1, j])
            label = None
            if edge_attr is not None:
                label = str(int(edge_attr[j].item()))
            g_nx.add_edge(u, v, label=label)
        graphs.append(Graph.from_networkx(g_nx))
    return graphs


DATASET_LOADERS: Dict[str, Callable[[], List[Graph]]] = {
    "toy": _load_toy,
    "mutag": _load_mutag,
    "enzymes": _load_enzymes,
}


def get_dataset(name: str, **kwargs) -> List[Graph]:
    """Load a dataset by name.

    Supported names include ``"toy"``, ``"mutag"`` and ``"enzymes"``.
    Names are case insensitive.

    Parameters
    ----------
    name: str
        Dataset identifier.
    **kwargs: dict
        Additional keyword arguments passed to the underlying loader. Not
        currently used.

    Returns
    -------
    List[Graph]
        List of graphs comprising the dataset.

    Raises
    ------
    KeyError
        If the dataset name is unknown.
    """
    key = name.lower()
    if key not in DATASET_LOADERS:
        raise KeyError(f"Unknown dataset '{name}'. Available datasets: {list(DATASET_LOADERS.keys())}")
    loader = DATASET_LOADERS[key]
    return loader()