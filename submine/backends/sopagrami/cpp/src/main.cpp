#include "alg.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <filesystem>
#include <functional>
#include <climits>

using namespace algo;
namespace fs = std::filesystem;



int main(int argc, char** argv){
    // Usage:
    //   run <graph.lg> [tau] [directed(0/1)] [sorted(0/1)] [threads]
    //
    // Defaults:
    //   tau=2, directed=0, sorted=1 (SoGraMi ordering), threads=4
   if (argc < 2){
    std::cerr
      << "Usage: run <graph.lg> [tau] [directed(0/1)] [sorted(0/1)] [threads]\n"
      << "               [dump_dir] [dump_images(0/1)] [max_images_per_vertex]\n"
      << "               [dump_emb(0/1)] [sample_limit]\n";
    return 1;
}


    const std::string path = argv[1];
    const int   tau      = (argc > 2 ? std::stoi(argv[2]) : 2);
    const bool  directed = (argc > 3 ? (std::stoi(argv[3]) != 0) : false);
    const bool  sorted   = (argc > 4 ? (std::stoi(argv[4]) != 0) : true);   // default: SoGraMi sorted
    const int   threads  = (argc > 5 ? std::stoi(argv[5]) : 4);             // default: 4

    DataGraph G;
    G.load_from_lg(path, directed);

    // Graph stats
    std::cout << "Graph loaded: |V|=" << G.vlabels.size() << ", |E|=";
    long long edge_count = 0;
    for (const auto& adj_list : G.adj) edge_count += (long long)adj_list.size();
    if (!directed) edge_count /= 2;
    std::cout << edge_count << "\n";

    // Params
    Params p;
    p.tau = tau;
    p.directed = directed;
    p.sorted_seeds = sorted;     // SoGraMi ordering toggle
    p.num_threads = threads;     // run_sopagrami  <=0 will default to all available
    p.compute_full_support = true;

    std::cout << "Settings: tau=" << p.tau
              << " directed=" << (p.directed?1:0)
              << " sorted=" << (p.sorted_seeds?1:0)
              << " threads=" << p.num_threads
              << "\n\n";

    // Run
    auto out = run_sopagrami(G, p);

    // Output
    std::cout << "Frequent patterns: " << out.frequent_patterns.size() << "\n";
    for (const auto& f : out.frequent_patterns){
        std::cout << "k=" << f.pat.vlab.size()
                  << " |E|=" << f.pat.pedges.size()
                  << " full=" << f.full_support
                  << " key=" << f.pat.key() << "\n";
    }
    //dump patterns to dir
   
    std::string dump_dir = (argc > 6 ? argv[6] : "");
    bool dump_images_csv = (argc > 7 ? (std::stoi(argv[7]) != 0) : false);
    int  max_images_per_vertex = (argc > 8 ? std::stoi(argv[8]) : 200);
    bool dump_sample_embeddings = (argc > 9 ? (std::stoi(argv[9]) != 0) : false);
    int  sample_limit = (argc > 10 ? std::stoi(argv[10]) : 50);

    if (!dump_dir.empty()){
        dump_patterns_to_dir(out, dump_dir, p.directed, G,
                            dump_images_csv, max_images_per_vertex,
                            dump_sample_embeddings, sample_limit);
        std::cout << "Wrote pattern files to: " << dump_dir
                << " (index.tsv, .lg, .dot"
                << (dump_images_csv ? ", .images.csv" : "")
                << (dump_sample_embeddings ? ", .emb.csv" : "")
                << ")\n";
    }

    return 0;
}
