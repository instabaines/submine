// python_bindings.cpp
#include "alg.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace algo;

// We’ll expose a very simple API: run on a .lg file and return patterns
// as a list of dict-like structures.
py::list run_on_lg_file(
    const std::string& path,
    int tau,
    bool directed,
    bool sorted_seeds,
    int num_threads,
    bool compute_full_support
) {
    DataGraph G;
    G.load_from_lg(path, directed);

    Params p;
    p.tau                = tau;
    p.directed           = directed;
    p.sorted_seeds       = sorted_seeds;
    p.num_threads        = num_threads;
    p.compute_full_support = compute_full_support;

    Output out = run_sopagrami(G, p);

    py::list py_patterns;
    for (const auto& f : out.frequent_patterns) {
        const Pattern& P = f.pat;

        py::dict d;
        d["node_labels"] = P.vlab;  // std::vector<std::string>

        py::list edges;
        for (const auto& e : P.pedges) {
            // (a, b, label, dir), dir: 0 undirected, 1 a->b
            edges.append(py::make_tuple(e.a, e.b, e.el, e.dir));
        }
        d["edges"]        = std::move(edges);
        d["full_support"] = f.full_support;
        d["key"]          = P.key();

        py_patterns.append(std::move(d));
    }
    return py_patterns;
}

PYBIND11_MODULE(sopagrami_cpp, m) {
    m.doc() = "pybind11 bindings for SoPaGraMi (C++17)";
    m.def(
        "run_on_lg_file",
        &run_on_lg_file,
        py::arg("path"),
        py::arg("tau")                 = 2,
        py::arg("directed")            = false,
        py::arg("sorted_seeds")        = true,
        py::arg("num_threads")         = 0,
        py::arg("compute_full_support")= true
    );
}
