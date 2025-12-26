#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <optional>
#include <tuple>
#include <limits>
#include <bitset>
#include <cstdint>
namespace algo {



// --- portable popcount64 ---
#if defined(_MSC_VER)
  #include <intrin.h>
  static inline int popcount64(unsigned long long x) {
  #if defined(_M_X64) || defined(_M_ARM64)
      return (int)__popcnt64(x);
  #else
      // 32-bit: split into two 32-bit halves
      return (int)(__popcnt((unsigned int)(x)) +
                   __popcnt((unsigned int)(x >> 32)));
  #endif
  }
#else
  static inline int popcount64(unsigned long long x) {
      return __builtin_popcountll(x);
  }
#endif


// ---------- DataGraph (your spec) ----------
struct Edge { int u, v; std::string label; };
struct Bitset {
    std::vector<uint64_t> w;
    void init(int n){ w.assign((n+63)>>6, 0ull); }
    inline void set(int i){ w[i>>6] |= (1ull<<(i&63)); }
    inline void reset(int i){ w[i>>6] &= ~(1ull<<(i&63)); }
    inline bool test(int i) const { return (w[i>>6] >> (i&63)) & 1ull; }
    inline bool any() const { for (auto x: w) if (x) return true; return false; }
    inline size_t count() const { size_t c=0; for (auto x: w) c += popcount64((unsigned long long)x); return c; }

    // this &= other
    inline void and_inplace(const Bitset& o){
        size_t m = w.size(); for (size_t i=0;i<m;++i) w[i] &= o.w[i];
    }
    // this &= ~other
    inline void andnot_inplace(const Bitset& o){
        size_t m = w.size(); for (size_t i=0;i<m;++i) w[i] &= ~o.w[i];
    }
    // (this & other).any() without allocating
    inline bool any_and(const Bitset& o) const {
        size_t m = w.size(); for (size_t i=0;i<m;++i) if (w[i] & o.w[i]) return true; return false;
    }

    // make a copy intersected with 'o'
    inline Bitset copy_and(const Bitset& o) const {
        Bitset t; t.w.resize(w.size());
        for (size_t i=0;i<w.size();++i) t.w[i] = w[i] & o.w[i];
        return t;
    }
};
struct DataGraph {
    bool directed = false;
    std::vector<std::string> vlabels; // node labels
    std::vector<std::vector<std::pair<int,std::string>>> adj, rev; // (nbr, label)
    std::vector<std::unordered_map<int, std::unordered_set<std::string>>> adj_set, rev_set; // fast has_edge
    std::unordered_map<std::string, std::unordered_set<int>> lab2nodes; // node-label -> nodes
    std::vector<std::unordered_map<std::string, std::vector<int>>> out_by_el, in_by_el;
    // Bitset indices (speedup for undirected and directed)
    std::unordered_map<std::string, Bitset> label_bits;                 // label -> nodes
    std::vector<std::unordered_map<std::string, Bitset>> adj_el_bits;   // per u: el -> bitset(neighbors via el)
        
    void load_from_lg(const std::string& path, bool as_directed);

    bool has_edge(int u, int v, const std::string& label) const {
        auto it = adj_set[u].find(v);
        if(it==adj_set[u].end()) return false;
        if(label.empty()) return !it->second.empty();
        return it->second.count(label)>0;
    }

    // Edge type key for frequent 1-edge seeds (SoGraMi)
    struct EdgeTypeKey {
        std::string lu, lv, el; int dirflag; // 0 undirected, 1 u->v
        bool operator==(EdgeTypeKey const& o) const {
            return lu==o.lu && lv==o.lv && el==o.el && dirflag==o.dirflag;
        }
        bool operator<(EdgeTypeKey const& o) const {
            if(lu!=o.lu) return lu<o.lu;
            if(lv!=o.lv) return lv<o.lv;
            if(el!=o.el) return el<o.el;
            return dirflag<o.dirflag;
        }
    };

    struct EdgeTypeStat { EdgeTypeKey key; int count; };

    // Python-parity: counts in FIRST-SEEN INSERTION ORDER while scanning adjacency
    std::vector<EdgeTypeStat> edge_type_counts_insertion_order() const;
};

// ---------- Algorithm API ----------
struct Params {
    int tau = 2;                 // frequency threshold
    bool directed = false;       // interpret graph as directed?
    bool sorted_seeds = true;    // SoGraMi ordering
    int num_threads = 0;         // 0 => use hardware_concurrency
    bool compute_full_support = true; // if false, use MNI only
};

struct Pattern {
    // Compact labeled pattern
    // nodes have labels; edges as (i,j,el,dirflag) with i<->j in [0..k-1]
    std::vector<std::string> vlab; // per pattern-node label
    struct PEdge { int a,b; std::string el; int dir; }; // dir: 0 undirected, 1 a->b
    std::vector<PEdge> pedges;

    // canonical key for de-duplication (label seq + sorted edges)
    std::string key() const;
    // in Pattern:
    std::vector<int> parent;   // parent[new_vertex] = anchor on RMP for forward edges, -1 for seed
    int rightmost;             // index of the last-added vertex

};

struct Found {
    Pattern pat;
    long long full_support = 0; // number of isomorphisms (or MNI if compute_full_support=false)
};

struct Output {
    std::vector<Found> frequent_patterns;
};
void dump_patterns_to_dir(
    const Output& out,
    const std::string& dump_dir,
    bool directed,
    const DataGraph& G,
    bool dump_images_csv,
    int  max_images_per_vertex,
    bool dump_sample_embeddings,
    int  sample_limit
);
Output run_sopagrami(const DataGraph& G, const Params& p);

} // namespace algo
