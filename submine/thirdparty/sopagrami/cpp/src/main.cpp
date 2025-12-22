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

// ---- utilities to reuse your edge-check logic ----
static inline bool ok_edge_map(const algo::DataGraph& G,
                               const algo::Pattern::PEdge& e,
                               int va, int vb,   // pattern endpoints
                               int ga, int gb)   // graph nodes mapped to (va,vb)
{
    if (e.dir == 1) {
        if (e.a == va && e.b == vb) return G.has_edge(ga, gb, e.el);
        if (e.a == vb && e.b == va) return G.has_edge(gb, ga, e.el);
        return true;
    } else {
        if ((e.a == va && e.b == vb) || (e.a == vb && e.b == va))
            return G.has_edge(ga, gb, e.el) || G.has_edge(gb, ga, e.el);
        return true;
    }
}

// Forward check: quick “does gi have any neighbor candidate for each incident edge?”
static bool forward_ok(const algo::DataGraph& G, const algo::Pattern& P,
                       int v, int gi,
                       const std::vector<int>& assign,
                       const std::vector<std::vector<int>>& dom)
{
    for (const auto& e : P.pedges){
        if (e.a!=v && e.b!=v) continue;
        int w = (e.a==v ? e.b : e.a);

        if (assign[w] != -1){
            if (e.a==v){ if (!ok_edge_map(G,e,e.a,e.b,gi,assign[w])) return false; }
            else       { if (!ok_edge_map(G,e,e.a,e.b,assign[w],gi)) return false; }
            continue;
        }

        bool okN = false;
        for (int gj : dom[w]){
            if (e.a==v){ if (ok_edge_map(G,e,e.a,e.b,gi,gj)) { okN=true; break; } }
            else       { if (ok_edge_map(G,e,e.a,e.b,gj,gi)) { okN=true; break; } }
        }
        if (!okN) return false;
    }
    return true;
}
// Try to find one *full embedding* with x_fixVar == fixNode.
// If found, returns true and fills 'assignment' (size k, graph node IDs).
static bool find_embedding_with_fixed(const algo::DataGraph& G,
                                      const algo::Pattern& P,
                                      int fixVar, int fixNode,
                                      std::vector<int>& assignment)
{
    const int k = (int)P.vlab.size();
    assignment.assign(k, -1);

    // Build label-consistent domains
    std::vector<std::vector<int>> dom(k);
    for (int i=0;i<k;++i){
        auto it = G.lab2nodes.find(P.vlab[i]);
        if (it == G.lab2nodes.end()) return false;
        dom[i].assign(it->second.begin(), it->second.end());
        if (dom[i].empty()) return false;
    }

    // Fix x_fixVar = fixNode
    assignment[fixVar] = fixNode;
    std::vector<char> used(G.vlabels.size(), 0);
    used[fixNode] = 1;

    auto choose = [&](){
        int best=-1, bestCnt=INT_MAX;
        for (int v=0; v<k; ++v){
            if (assignment[v]!=-1) continue;
            int cnt=0;
            for (int gi : dom[v]){
                if (used[gi]) continue;
                if (forward_ok(G,P,v,gi,assignment,dom)){ ++cnt; if (cnt>=bestCnt) break; }
            }
            if (cnt < bestCnt){ best=v; bestCnt=cnt; }
        }
        return best;
    };

    std::function<bool()> dfs = [&](){
        for (int i=0;i<k;++i) if (assignment[i]==-1) goto not_done;
        return true;
      not_done:
        int v = choose(); if (v==-1) return false;
        for (int gi : dom[v]){
            if (used[gi]) continue;
            if (!forward_ok(G,P,v,gi,assignment,dom)) continue;
            assignment[v]=gi; used[gi]=1;
            if (dfs()) return true;
            used[gi]=0; assignment[v]=-1;
        }
        return false;
    };

    return dfs();
}


// For each pattern vertex i, collect up to `max_per_vertex` graph node IDs
// that participate in at least one full embedding (MNI “image set”).
// If max_per_vertex < 0 => no cap.
static std::vector<std::vector<int>>
collect_mni_image_sets(const algo::DataGraph& G,
                       const algo::Pattern& P,
                       int max_per_vertex = 100)
{
    const int k = (int)P.vlab.size();
    std::vector<std::vector<int>> images(k);

    // Domains by label
    std::vector<std::vector<int>> dom(k);
    for (int i=0;i<k;++i){
        auto it = G.lab2nodes.find(P.vlab[i]);
        if (it == G.lab2nodes.end()) return images;
        dom[i].assign(it->second.begin(), it->second.end());
    }

    // For each pattern variable v, test each u in dom[v] by trying to find one embedding
    for (int v=0; v<k; ++v){
        int kept = 0;
        for (int u : dom[v]){
            std::vector<int> a;
            if (find_embedding_with_fixed(G, P, v, u, a)){
                images[v].push_back(u);
                ++kept;
                if (max_per_vertex >= 0 && kept >= max_per_vertex) break;
            }
        }
    }
    return images;
}



// --- NEW: per-vertex images CSV (patternIndex, graphNodeId) ---
static void write_pattern_images_csv(const algo::Pattern& P,
                                     const std::vector<std::vector<int>>& images,
                                     const std::string& path_csv)
{
    std::ofstream out(path_csv);
    if (!out) return;
    out << "pattern_vertex,graph_node_id\n";
    for (size_t i=0;i<images.size();++i){
        for (int u : images[i]){
            out << i << "," << u << "\n";
        }
    }
}


// --- NEW: sample embeddings CSV (one row per embedding, columns are pattern vertex order) ---
static void write_sample_embeddings_csv(const algo::Pattern& P,
                                        const std::vector<std::vector<int>>& emb,
                                        const std::string& path_csv)
{
    std::ofstream out(path_csv);
    if (!out) return;
    // header
    out << "emb_id";
    for (size_t i=0;i<P.vlab.size();++i) out << ",v" << i;
    out << "\n";
    for (size_t i=0;i<emb.size();++i){
        out << i;
        for (int id : emb[i]) out << "," << id;
        out << "\n";
    }
}





// _____________________________________________________
static std::string sanitize_dot(const std::string& s){
    std::string t; t.reserve(s.size()*2);
    for (char c: s){
        if (c=='"' || c=='\\') t.push_back('\\');
        t.push_back(c);
    }
    return t;
}

static void write_pattern_as_lg(const algo::Pattern& P, const std::string& path){
    std::ofstream out(path);
    if (!out) return;
    for (size_t i=0;i<P.vlab.size();++i) out << "v " << i << " " << P.vlab[i] << "\n";
    for (const auto& e : P.pedges)       out << "e " << e.a << " " << e.b << " " << e.el << "\n";
}

static void write_pattern_as_dot(const algo::Pattern& P, bool directed, const std::string& path){
    std::ofstream out(path);
    if (!out) return;
    out << (directed ? "digraph G {\n" : "graph G {\n");
    // nodes
    for (size_t i=0;i<P.vlab.size();++i){
        out << "  " << i << " [shape=circle,label=\"" << sanitize_dot(P.vlab[i]) << "\"];\n";
    }
    // edges
    for (const auto& e : P.pedges){
        const bool use_arrow = directed || e.dir==1;
        out << "  " << e.a << (use_arrow ? " -> " : " -- ") << e.b
            << " [label=\"" << sanitize_dot(e.el) << "\"];\n";
    }
    out << "}\n";
}

// ---- Dump for static miner output ----
static void dump_patterns_to_dir(const algo::Output& out,
                                 const std::string& dir,
                                 bool directed,
                                 const algo::DataGraph& G,
                                 bool dump_images_csv = false,
                                 int  max_images_per_vertex = 200,
                                 bool dump_sample_embeddings = false,
                                 int  sample_limit = 50)
{
    namespace fs = std::filesystem;
    fs::create_directories(dir);
    std::ofstream idx(dir + "/index.tsv");
    idx << "id\tk\tm\tfull_support\tkey\tlg_path\tdot_path\n";

    for (size_t i=0; i<out.frequent_patterns.size(); ++i){
        const auto& f = out.frequent_patterns[i];
        const size_t k = f.pat.vlab.size();
        const size_t m = f.pat.pedges.size();

        std::string base = dir + "/pat_" + std::to_string(i)
                         + "_k" + std::to_string(k)
                         + "_e" + std::to_string(m)
                         + "_full" + std::to_string(f.full_support);
        std::string lgp  = base + ".lg";
        std::string dotp = base + ".dot";

        // always write shape artifacts
        write_pattern_as_lg (f.pat, lgp);
        write_pattern_as_dot(f.pat, directed, dotp);

        // optionally: image sets (can be heavy)
        if (dump_images_csv){
            auto images = collect_mni_image_sets(G, f.pat, max_images_per_vertex);
            write_pattern_images_csv(f.pat, images, base + ".images.csv");
        }

        // optionally: sample full embeddings (disabled in your current code; left stub)
        if (dump_sample_embeddings){
            std::vector<std::vector<int>> samples;
            // enumerate_embeddings(G, f.pat, sample_limit, samples); // if/when you implement
            write_sample_embeddings_csv(f.pat, samples, base + ".emb.csv");
        }

        idx << i << '\t' << k << '\t' << m << '\t'
            << f.full_support << '\t' << f.pat.key()
            << '\t' << lgp << '\t' << dotp << "\n";
    }
}



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
