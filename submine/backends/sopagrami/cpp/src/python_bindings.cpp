// python_bindings.cpp
#include "alg.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace algo;

// Expose: run on a .lg file, optionally dump pattern artifacts to a directory,
// and return patterns as a list of dicts.
py::list run_on_lg_file(
    const std::string& path,
    int  tau,
    bool directed,
    bool sorted_seeds,
    int  num_threads,
    bool compute_full_support,
    const std::string& out_dir,
    bool dump_images_csv,
    int  max_images_per_vertex,
    bool dump_sample_embeddings,
    int  sample_limit
) {
    DataGraph G;
    G.load_from_lg(path, directed);

    Params p;
    p.tau                  = tau;
    p.directed             = directed;
    p.sorted_seeds         = sorted_seeds;
    p.num_threads          = num_threads;
    p.compute_full_support = compute_full_support;

    Output out = run_sopagrami(G, p);

    // Optional side-effect: dump pattern files to directory
    if (!out_dir.empty()) {
        dump_patterns_to_dir(
            out,
            out_dir,
            p.directed,
            G,
            dump_images_csv,
            max_images_per_vertex,
            dump_sample_embeddings,
            sample_limit
        );
    }

    // Return patterns to Python
    py::list py_patterns;

    for (const auto& f : out.frequent_patterns) {
        const Pattern& P = f.pat;

        py::dict d;
        d["node_labels"] = P.vlab; // std::vector<std::string>

        py::list edges;
        for (const auto& e : P.pedges) {
            // (a, b, label, dir) dir: 0 undirected, 1 a->b (per your comment)
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
        py::arg("tau")                  = 2,
        py::arg("directed")             = false,
        py::arg("sorted_seeds")         = true,
        py::arg("num_threads")          = 0,
        py::arg("compute_full_support") = true,

        // dump-related args 
        py::arg("out_dir")              = std::string("result"),
        py::arg("dump_images_csv")      = false,
        py::arg("max_images_per_vertex")= 200,
        py::arg("dump_sample_embeddings")= false,
        py::arg("sample_limit")         = 50,

        R"doc(
Run SoPaGraMi on an input .lg graph.

Parameters
----------
path : str
    Path to input .lg file.
tau : int, default=2
directed : bool, default=False
sorted_seeds : bool, default=True
num_threads : int, default=0
    0 means "use default / auto" as implemented in C++ core.
compute_full_support : bool, default=True

out_dir : str, default=""
    If non-empty, dumps pattern artifacts to this directory:
    index.tsv, per-pattern .lg, .dot, plus optional .images.csv and .emb.csv.
dump_images_csv : bool, default=False
max_images_per_vertex : int, default=200
dump_sample_embeddings : bool, default=False
sample_limit : int, default=50

Returns
-------
list[dict]
    Each dict contains: node_labels, edges, full_support, key.
)doc"
    );
}
