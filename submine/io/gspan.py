# submine/io/gspan.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, Hashable, Iterable, List, Tuple

from ..core.graph import Graph
from ..errors import SubmineInputError
from ..utils.checks import iter_text_lines



from pathlib import Path
from typing import Iterable, List, Optional



def read_gspan_dataset(path: Path | str) -> List[Graph]:
    """
    Read a gSpan-formatted dataset file and return a list of Graph objects.

    Expected format:

        t # N         # start of N-th graph
        v M L         # vertex M has label L
        e P Q L       # edge (P, Q) has label L
        ...
        t # -1        # end of file sentinel (required by some gSpan impls)

    Notes
    -----
    - Vertex ids within each graph are assumed to be integers (0..n-1).
    - Labels are read as integers (you can remap them later if desired).
    - Edges are treated as undirected; we store them as (min(u, v), max(u, v))
      and avoid duplicates.
    """
    path = Path(path)

    graphs: List[Graph] = []

    current_nodes: Optional[list[int]] = None
    current_node_labels: Optional[dict[int, int]] = None
    current_edges: Optional[list[tuple[int, int]]] = None
    current_edge_labels: Optional[dict[tuple[int, int], int]] = None

    def flush_current_graph():
        nonlocal current_nodes, current_node_labels, current_edges, current_edge_labels
        if current_nodes is None:
            return
        graphs.append(
            Graph(
                nodes=current_nodes,
                edges=current_edges,
                node_labels=current_node_labels,
                edge_labels=current_edge_labels,
            )
        )
        current_nodes = None
        current_node_labels = None
        current_edges = None
        current_edge_labels = None

    for raw_line in iter_text_lines(path):
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split()

        rec_type = parts[0]

        # Graph header: t # gid
        if rec_type == "t":
            # flush previous graph if any
            if len(parts) >= 3:
                gid = int(parts[2])
            else:
                raise ValueError(f"Malformed 't' line: {line!r}")

            if gid == -1:
                # End-of-dataset sentinel: flush the last graph and stop
                flush_current_graph()
                break

            # start a new graph
            flush_current_graph()
            current_nodes = []
            current_node_labels = {}
            current_edges = []
            current_edge_labels = {}

        elif rec_type == "v":
            if current_nodes is None:
                raise ValueError(f"Vertex line outside of any graph: {line!r}")
            if len(parts) < 3:
                raise ValueError(f"Malformed 'v' line: {line!r}")
            vid = int(parts[1])
            lbl = int(parts[2])
            current_nodes.append(vid)
            current_node_labels[vid] = lbl  # type: ignore[arg-type]

        elif rec_type == "e":
            if current_edges is None:
                raise ValueError(f"Edge line outside of any graph: {line!r}")
            if len(parts) < 4:
                raise ValueError(f"Malformed 'e' line: {line!r}")
            u = int(parts[1])
            v = int(parts[2])
            lbl = int(parts[3])

            # treat as undirected, avoid duplicates
            if u <= v:
                key = (u, v)
            else:
                key = (v, u)

            if key not in current_edge_labels:  # type: ignore[operator]
                current_edges.append(key)       # type: ignore[arg-type]
                current_edge_labels[key] = lbl  # type: ignore[index]

        else:
            # Unknown record; you can choose to ignore or raise
            raise SubmineInputError(f"Unknown record type '{rec_type}' in line: {line!r}")

    # If file ended without the required terminator, still return what we have.
    flush_current_graph()
    return graphs

def convert_gspan_graph(gspan_g) -> Graph:
    """
    Convert a vendored gSpan graph.Graph object into submine.core.graph.Graph.

    Assumes:
      - gspan_g.vertices is a dict {vid: Vertex}
      - Vertex has: vid, vlb, edges (dict[to_vid, Edge])
      - Edge has: eid, frm, to, elb
    """
    nodes = []
    node_labels = {}
    edges = []
    edge_labels = {}

    # 1. Nodes + labels
    for vid, v in gspan_g.vertices.items():
        nodes.append(vid)
        node_labels[vid] = v.vlb

    # 2. Edges (avoid duplicates in undirected graphs)
    seen = set()
    for vid, v in gspan_g.vertices.items():
        # v.edges is a dict: {to_vid: Edge}
        for to, e in v.edges.items():
            u, w = e.frm, e.to
            # canonicalize for undirected graph
            key = (u, w) if u <= w else (w, u)
            if key in seen:
                continue
            seen.add(key)
            edges.append(key)
            edge_labels[key] = e.elb

    return Graph(
        nodes=nodes,
        edges=edges,
        node_labels=node_labels,
        edge_labels=edge_labels,
    )


def _build_label_maps(graphs: List[Graph]):
    """
    Map arbitrary node/edge labels to consecutive ints >= 2,
    because gSpan forbids 0 and 1.
    """
    node_label_map: Dict[Hashable, int] = {}
    edge_label_map: Dict[Hashable, int] = {}

    next_node_label = 2
    next_edge_label = 2

    for G in graphs:
        # Node labels
        if G.node_labels is not None:
            for nid in G.nodes:
                lbl = G.node_labels.get(nid, None)
                if lbl is None:
                    continue
                if lbl not in node_label_map:
                    node_label_map[lbl] = next_node_label
                    next_node_label += 1

        # Edge labels
        if G.edge_labels is not None:
            for e, lbl in G.edge_labels.items():
                if lbl not in edge_label_map:
                    edge_label_map[lbl] = next_edge_label
                    next_edge_label += 1

    # Fallback: if there are unlabeled nodes/edges, give them a default label
    if not node_label_map:
        node_label_map["__default_node__"] = 2
    if not edge_label_map:
        edge_label_map["__default_edge__"] = 2

    return node_label_map, edge_label_map


def write_gspan_dataset(graphs: Iterable[Graph], path: Path) -> None:
    """
    Write a list of Graph objects to a gSpan-compatible file.

    Format:
        t # N           -> N-th graph
        v M L           -> vertex M has label L
        e P Q L         -> edge (P, Q) has label L
        ...
        t # -1          -> end of file

    NOTE:
      - Vertex ids must be 0..n-1 *within each graph*.
      - All labels must be integers >= 2 (we map them if needed).
    """
    
    graphs = list(graphs)
    node_label_map, edge_label_map = _build_label_maps(graphs)

    with path.open("w") as f:
        for gid, G in enumerate(graphs):
            f.write(f"t # {gid}\n")

            # Remap node ids to 0..n-1 locally
            id_map = {orig_id: new_id for new_id, orig_id in enumerate(G.nodes)}

            # Write vertices
            for orig_id in G.nodes:
                new_id = id_map[orig_id]
                if G.node_labels is None:
                    # use default label
                    label_key = "__default_node__"
                else:
                    raw_lbl = G.node_labels.get(orig_id, "__default_node__")
                    label_key = raw_lbl if raw_lbl in node_label_map else "__default_node__"

                lbl_int = node_label_map[label_key]
                f.write(f"v {new_id} {lbl_int}\n")

            # Write edges
            for e in G.edges:
                if len(e) == 2:
                    u, v = e
                    raw_elbl = "__default_edge__"
                elif len(e) == 3:
                    u, v, raw_elbl = e
                else:
                    raise ValueError(f"Edge tuple must be (u,v) or (u,v,label), got {e!r}")

                u_new = id_map[u]
                v_new = id_map[v]

                if G.edge_labels is not None:
                    raw_elbl = G.edge_labels.get((u, v), raw_elbl)

                label_key = raw_elbl if raw_elbl in edge_label_map else "__default_edge__"
                elbl_int = edge_label_map[label_key]

                f.write(f"e {u_new} {v_new} {elbl_int}\n")

        # Required terminator for this implementation
        f.write("t # -1\n")