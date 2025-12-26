#include "alg.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <functional>
#include <climits>
namespace fs = std::filesystem;

namespace algo {
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



// per-vertex images CSV (patternIndex, graphNodeId) ---
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

void dump_patterns_to_dir(
    const Output& out,
    const std::string& dump_dir,
    bool directed,
    const DataGraph& G,
    bool dump_images_csv,
    int  max_images_per_vertex,
    bool dump_sample_embeddings,
    int  sample_limit
) {
    fs::create_directories(dump_dir);

    // ---- index.tsv ----
    std::ofstream idx(fs::path(dump_dir) / "index.tsv");
    idx << "id\tk\tm\tfull_support\tkey\tlg_path\tdot_path\n";

    for (size_t i=0; i<out.frequent_patterns.size(); ++i){
        const auto& f = out.frequent_patterns[i];
        const size_t k = f.pat.vlab.size();
        const size_t m = f.pat.pedges.size();

        std::string base = dump_dir + "/pat_" + std::to_string(i)
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
            // enumerate_embeddings(G, f.pat, sample_limit, samples); // not implelemnted yet
            write_sample_embeddings_csv(f.pat, samples, base + ".emb.csv");
        }

        idx << i << '\t' << k << '\t' << m << '\t'
            << f.full_support << '\t' << f.pat.key()
            << '\t' << lgp << '\t' << dotp << "\n";
    }
}

} // namespace algo
