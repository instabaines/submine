#include <./pybind11/pybind11.h>
#include <./pybind11/stl.h>

#include <sstream>
#include <unordered_set>
#include <vector>
#include <utility>
#include <cstdint>
#include <algorithm>

#include "gspan.h"

namespace py = pybind11;

// Hash for pair<int,int>
struct PairHash {
    std::size_t operator()(const std::pair<int,int>& p) const noexcept {
        return (static_cast<std::size_t>(p.first) << 32) ^ static_cast<std::size_t>(p.second);
    }
};

static py::dict graph_to_dict(const GSPAN::Graph& g,
                              unsigned int support,
                              const GSPAN::Projected* projected,
                              bool directed)
{
    py::dict d;
    const int n = static_cast<int>(g.size());

    std::vector<int> nodes;
    nodes.reserve(n);

    std::vector<int> node_labels;
    node_labels.reserve(n);

    for (int i = 0; i < n; ++i) {
        nodes.push_back(i);
        node_labels.push_back(static_cast<int>(g[i].label));
    }

    std::vector<std::pair<int,int>> edges;
    std::vector<int> edge_labels;

    // IMPORTANT: raw extraction, no dedup
    for (int u = 0; u < n; ++u) {
        for (const auto& e : g[u].edge) {
            edges.emplace_back(u, static_cast<int>(e.to));
            edge_labels.push_back(static_cast<int>(e.elabel));
        }
    }


    d["nodes"] = std::move(nodes);
    d["edges"] = std::move(edges);
    d["node_labels"] = std::move(node_labels);
    d["edge_labels"] = std::move(edge_labels);
    d["support"] = support;

    // Graph IDs (gid) where pattern occurs.
    // In your code: PDFS.id is "ID of the original input graph"
    if (projected) {
        std::unordered_set<unsigned int> uniq;
        uniq.reserve(projected->size());

        std::vector<unsigned int> gids;
        gids.reserve(projected->size());

        for (const auto& p : *projected) {
            unsigned int gid = p.id;   // <-- this is your gid
            if (uniq.insert(gid).second) gids.push_back(gid);
        }

        std::sort(gids.begin(), gids.end());
        d["graph_ids"] = std::move(gids);
    } else {
        d["graph_ids"] = py::none();
    }

    return d;
}


static std::vector<py::dict> mine_from_string(const std::string& gspan_data,
                                              unsigned int minsup,
                                              unsigned int maxpat_min,
                                              unsigned int maxpat_max,
                                              bool enc,
                                              bool where,
                                              bool directed)
{
    std::istringstream is(gspan_data);
    std::ostringstream null_out; // discard textual output by default

    GSPAN::gSpan miner;
    std::vector<py::dict> results;
    results.reserve(1024);

    miner.set_callback([&](const GSPAN::Graph& pattern,
                           unsigned int sup,
                           const GSPAN::Projected* projected)
    {
       results.push_back(graph_to_dict(pattern, sup, projected, directed));

    });

    miner.run(is, null_out, minsup, maxpat_min, maxpat_max, enc, where, directed);
    return results;
}

PYBIND11_MODULE(gspan_cpp, m) {
    m.doc() = "Pure C++ gSpan bindings (structured results)";

    m.def("mine_from_string",
          &mine_from_string,
          py::arg("gspan_data"),
          py::arg("minsup") = 1,
          py::arg("maxpat_min") = 0,
          py::arg("maxpat_max") = 0xffffffffu,
          py::arg("enc") = false,
          py::arg("where") = false,
          py::arg("directed") = false,
          R"pbdoc(
Mine frequent subgraphs from a gSpan text dataset passed as a string.

Returns: list of dicts:
  - nodes: [0..n-1]
  - edges: list of (u,v)
  - node_labels: list aligned to nodes
  - edge_labels: list aligned to edges
  - support: int
  - graph_ids: list of graph IDs where pattern occurs, or None
)pbdoc");
}
