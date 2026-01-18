// alg.cpp
#include "alg.hpp"

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <climits>
#include <map>
#include <tuple>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace algo {

// =========================== DataGraph ===============================

static inline string encKey(const DataGraph::EdgeTypeKey& k){
    return k.lu + '\t' + k.lv + '\t' + k.el + '\t' + char('0' + k.dirflag);
}

// ---------------------------------------------------------------------
// 1. Build Indices: Create Integer IDs and Bitsets for O(1) Access
// ---------------------------------------------------------------------
void DataGraph::build_indices() {
    int n = (int)vlabels.size();

    // A. Map String Labels to Integer IDs
    label_to_id.clear();
    id_to_label.clear();
    int next_id = 0;

    auto get_id = [&](const std::string& s) {
        if(label_to_id.find(s) == label_to_id.end()) {
            label_to_id[s] = next_id++;
            id_to_label.push_back(s);
        }
        return label_to_id[s];
    };

    // Scan graph to build ID map
    for(const auto& list : adj) {
        for(const auto& e : list) get_id(e.second);
    }

    // B. Initialize Vector-based Bitsets [Num_Labels][Num_Nodes]
    int num_labels = next_id;
    adj_bits.assign(num_labels, std::vector<Bitset>(n));
    rev_bits.assign(num_labels, std::vector<Bitset>(n));

    for(int l=0; l<num_labels; ++l) {
        for(int u=0; u<n; ++u) {
            adj_bits[l][u].init(n);
            rev_bits[l][u].init(n);
        }
    }

    // C. Populate Bitsets & Lab2Nodes
    lab2nodes.clear();
    for(int i=0; i<n; ++i) lab2nodes[vlabels[i]].insert(i);

    for(int u=0; u<n; ++u) {
        for(const auto& e : adj[u]) {
            int v = e.first;
            int l_id = label_to_id[e.second];

            adj_bits[l_id][u].set(v);
            rev_bits[l_id][v].set(u);
        }
    }
}

void DataGraph::load_from_lg(const std::string& path, bool as_directed){
    directed = as_directed;
    std::ifstream fin(path);
    if(!fin){ std::cerr<<"Cannot open "<<path<<"\n"; std::exit(1); }

    struct EdgeRec { int u,v; std::string el; };
    std::unordered_map<int,std::string> vlab_map;
    std::vector<EdgeRec> edges;

    std::string line;
    while(std::getline(fin,line)){
        if(line.empty()) continue;
        std::stringstream ss(line);
        std::string tag; ss>>tag;

        if(tag=="v"||tag=="V"){
            int idx; std::string lab; ss>>idx>>lab;
            vlab_map[idx]=lab;
        }else if(tag=="e"||tag=="E"){
            int u,v; std::string elab; ss>>u>>v;
            if(!(ss>>elab)) elab = "";
            if(!elab.empty()){
                bool numlike=true;
                for(unsigned char c : elab){
                    if(!(std::isdigit(c)||c=='.'||c=='-'||c=='+')){ numlike=false; break; }
                }
                if(numlike){
                    try{ double d = std::stod(elab); long long iv = (long long)d; elab = std::to_string(iv); }catch(...){}
                }
            }
            edges.push_back({u,v,elab});
        }
    }

    const int n = vlab_map.empty()? 0 :
        (std::max_element(vlab_map.begin(), vlab_map.end(),
            [](auto&a, auto&b){ return a.first < b.first; })->first + 1);

    vlabels.assign(n,"");
    for(int i=0;i<n;++i) if(vlab_map.count(i)) vlabels[i]=vlab_map[i];

    adj.assign(n,{}); rev.assign(n,{});
    adj_set.assign(n,{}); rev_set.assign(n,{});

    for(const auto& e: edges){
        adj[e.u].push_back({e.v,e.el});
        adj_set[e.u][e.v].insert(e.el);
        rev[e.v].push_back({e.u,e.el});
        rev_set[e.v][e.u].insert(e.el);
        if(!directed){
            adj[e.v].push_back({e.u,e.el});
            adj_set[e.v][e.u].insert(e.el);
            rev[e.u].push_back({e.v,e.el});
            rev_set[e.u][e.v].insert(e.el);
        }
    }
    // Call the builder to setup optimized structures
    build_indices();
}

// ============================ Pattern =================================

string Pattern::key() const {
    string s; s.reserve(vlab.size()*8 + pedges.size()*16);
    s += "V:";
    for (size_t i=0;i<vlab.size();++i){ s += vlab[i]; s += '|'; }

    vector<tuple<int,int,int,string>> es; es.reserve(pedges.size());
    for(const auto& e: pedges){
        int a=e.a,b=e.b;
        int dcode = 0; // 0 undirected, 1 a->b oriented to min(a,b), 2 b->a oriented
        if(e.dir==1) dcode = (a<b ? 1 : 2);
        es.emplace_back(min(a,b), max(a,b), dcode, e.el);
    }
    sort(es.begin(), es.end());
    s += "E:";
    for(auto& t: es){
        s += to_string(get<0>(t))+"-"+to_string(get<1>(t))+"-"+to_string(get<2>(t))+"-"+get<3>(t)+"|";
    }
    return s;
}


static std::string canonical_key(const Pattern& S){
    const int k = (int)S.vlab.size();
    if (k <= 1) return S.key();

    // 1. Compute Vertex Signatures (Label + Degree + Sorted Neighbor Labels)
    // This distinguishes nodes that look identical (same label) but are structurally different.
    std::vector<std::pair<std::string, int>> node_sigs(k);

    // Pre-build adjacency for the PATTERN (not the data graph)
    std::vector<std::vector<std::string>> pat_adj(k);
    std::vector<int> degrees(k, 0);

    for(const auto& e : S.pedges) {
        degrees[e.a]++;
        degrees[e.b]++;
        pat_adj[e.a].push_back(S.vlab[e.b] + e.el); // Neighbor label + Edge label
        pat_adj[e.b].push_back(S.vlab[e.a] + e.el);
    }

    for(int i=0; i<k; ++i) {
        std::sort(pat_adj[i].begin(), pat_adj[i].end());
        std::string neighbor_str;
        for(const auto& s : pat_adj[i]) neighbor_str += s;

        // Signature: "Label_Degree_NeighborString"
        // We use a pair {Signature, OriginalIndex} to sort/group
        node_sigs[i] = { S.vlab[i] + "_" + std::to_string(degrees[i]) + "_" + neighbor_str, i };
    }

    // 2. Group by Unique Signature (Refinement)
    // Instead of grouping just by "Label", we group by the detailed signature.
    std::map<std::string, std::vector<int>> groups;
    for(int i=0; i<k; ++i) {
        groups[node_sigs[i].first].push_back(i);
    }

    // Sort the original indices within groups to ensure determinism
    for(auto& kv : groups) {
        std::sort(kv.second.begin(), kv.second.end());
    }

    // 3. Deterministic Iteration Order
    std::vector<std::string> sorted_sig_keys;
    for(const auto& kv : groups) sorted_sig_keys.push_back(kv.first);

    // 4. Backtracking Permutation (Now much faster)
    auto encode_with_order = [&](const std::vector<int>& order)->std::string{
        std::vector<int> pos(k);
        for (int i=0;i<k;++i) pos[order[i]] = i;

        std::string s; s.reserve(k*8 + S.pedges.size()*16);
        s += "V:";
        for (int i=0;i<k;++i){ s += S.vlab[order[i]]; s += '|'; }

        std::vector<std::tuple<int,int,int,std::string>> es;
        es.reserve(S.pedges.size());
        for (const auto& e : S.pedges){
            int a = pos[e.a], b = pos[e.b];
            int dcode = 0;
            if (e.dir==1) dcode = (a<b ? 1 : 2);
            es.emplace_back(std::min(a,b), std::max(a,b), dcode, e.el);
        }
        std::sort(es.begin(), es.end());

        s += "E:";
        for (auto& t : es){
            s += std::to_string(std::get<0>(t)) + "-" +
                 std::to_string(std::get<1>(t)) + "-" +
                 std::to_string(std::get<2>(t)) + "-" +
                 std::get<3>(t) + "|";
        }
        return s;
    };

    std::string best; bool have=false;
    std::vector<int> current; current.reserve(k);

    std::function<void(int)> dfs = [&](int group_idx){
        if (group_idx == (int)sorted_sig_keys.size()){
            std::string code = encode_with_order(current);
            if (!have || code < best){ best = std::move(code); have = true; }
            return;
        }

        // Permute only within the refined group
        const auto& g = groups[sorted_sig_keys[group_idx]];
        std::vector<int> perm = g;

        do {
            size_t old_size = current.size();
            current.insert(current.end(), perm.begin(), perm.end());

            // Optimization: If partial sequence is already "worse" than best, we could prune.
            // But 'encode_with_order' needs full sequence.
            // For now, refinement is the main speedup.

            dfs(group_idx + 1);

            current.resize(old_size);
        } while (std::next_permutation(perm.begin(), perm.end()));
    };

    dfs(0);
    return best;
}

// ======================= Seeds (1-edge) ================================

struct SeedInfo {
    DataGraph::EdgeTypeKey key;
    int mni;           // MNI support for the 1-edge pattern
    long long full;    // number of edges of that type
};

static vector<SeedInfo> compute_frequent_edge_seeds(const DataGraph& G, int tau){
    using K = DataGraph::EdgeTypeKey;
    struct AccDir { unordered_set<int> L, R; long long full=0; K key; };
    struct AccEq  { unordered_set<int> S;   long long full=0; K key; };

    unordered_map<string, AccDir> acc_lr; acc_lr.reserve(1<<14);
    unordered_map<string, AccEq>  acc_eq; acc_eq.reserve(1<<14);
    const int n = (int)G.vlabels.size();

    if (G.directed){
        for (int u=0; u<n; ++u){
            const string& lu = G.vlabels[u];
            for (auto [v, el] : G.adj[u]){
                const string& lv = G.vlabels[v];
                K k{lu, lv, el, 1};
                auto &A = acc_lr[ encKey(k) ];
                if (A.full == 0) A.key = k;
                A.L.insert(u); A.R.insert(v); A.full += 1;
            }
        }
    } else {
        for (int u=0; u<n; ++u){
            const string& lu = G.vlabels[u];
            for (auto [v, el] : G.adj[u]){
                if (u > v) continue;
                const string& lv = G.vlabels[v];
                if (lu == lv){
                    K k{lu, lv, el, 0};
                    auto &E = acc_eq[ encKey(k) ];
                    if (E.full == 0) E.key = k;
                    E.S.insert(u); E.S.insert(v); E.full += 1;
                } else {
                    K k;
                    int leftNode, rightNode;
                    if (lu <= lv) { k = {lu, lv, el, 0}; leftNode=u; rightNode=v; }
                    else          { k = {lv, lu, el, 0}; leftNode=v; rightNode=u; }
                    auto &A = acc_lr[ encKey(k) ];
                    if (A.full == 0) A.key = k;
                    A.L.insert(leftNode); A.R.insert(rightNode); A.full += 1;
                }
            }
        }
    }

    vector<SeedInfo> seeds;
    for (auto &kv : acc_lr){
        auto &A = kv.second;
        int mni = std::min((int)A.L.size(), (int)A.R.size());
        if (mni >= tau) seeds.push_back({A.key, mni, A.full});
    }
    for (auto &kv : acc_eq){
        auto &E = kv.second;
        int mni = (int)E.S.size();
        if (mni >= tau) seeds.push_back({E.key, mni, E.full});
    }
    return seeds;
}

static unordered_map<string,int> build_seed_mni_map(const vector<SeedInfo>& seeds){
    unordered_map<string,int> m; m.reserve(seeds.size()*2);
    for (auto &s : seeds) m[encKey(s.key)] = s.mni;
    return m;
}

// ==================== MNI & Support Calculation ===========================

// ---------------------------------------------------------------------
// 2. Helper: Fast O(1) Edge Check using Bitsets
// ---------------------------------------------------------------------
static inline bool check_edge_fast(const DataGraph& G, int dir, int label_id, int u, int v) {
    if (label_id < 0 || label_id >= (int)G.adj_bits.size()) return false;
    // G.adj_bits is symmetric for undirected graphs, so this logic covers both
    return G.adj_bits[label_id][u].test(v);
}

// ---------------------------------------------------------------------
// 3. Filter Domains (AC-3): Use Bitsets for Speed
// ---------------------------------------------------------------------
static void filter_domains_by_local_constraints(const DataGraph& G,
                                                const Pattern& P,
                                                vector<vector<int>>& dom_vecs)
{
    const int n = (int)G.vlabels.size();
    const int k = (int)P.vlab.size();

    // Convert vector domains to Bitsets
    vector<Bitset> dom(k);
    for(int i=0; i<k; ++i){
        dom[i].init(n);
        for(int u : dom_vecs[i]) dom[i].set(u);
    }

    bool changed = true;
    while(changed) {
        changed = false;
        for(const auto& e : P.pedges) {
            int u = e.a;
            int v = e.b;

            // Optimization: Look up Integer ID for edge label
            auto it_lab = G.label_to_id.find(e.el);
            int label_id = (it_lab != G.label_to_id.end()) ? it_lab->second : -1;

            if (label_id == -1) {
                 dom[u].init(n); dom[v].init(n); goto convert_back;
            }

            // Filter U based on V
            for (int cand_u = 0; cand_u < n; ++cand_u) {
                if (!dom[u].test(cand_u)) continue;

                bool keep = false;
                // AC-3 Logic: Does cand_u have ANY neighbor in dom[v] via label_id?
                if(G.adj_bits[label_id][cand_u].any_and(dom[v])) keep = true;

                if (!keep) { dom[u].reset(cand_u); changed = true; }
            }
            if (dom[u].count() == 0) goto convert_back;

            // Filter V based on U
            for (int cand_v = 0; cand_v < n; ++cand_v) {
                if (!dom[v].test(cand_v)) continue;

                bool keep = false;
                if (G.directed && e.dir == 1) {
                    // Check Incoming (rev_bits)
                    if (G.rev_bits[label_id][cand_v].any_and(dom[u])) keep = true;
                } else {
                    if (G.adj_bits[label_id][cand_v].any_and(dom[u])) keep = true;
                }

                if (!keep) { dom[v].reset(cand_v); changed = true; }
            }
            if (dom[v].count() == 0) goto convert_back;
        }
    }

convert_back:
    for(int i=0; i<k; ++i){
        dom_vecs[i].clear();
        for(int u=0; u<n; ++u){
            if(dom[i].test(u)) dom_vecs[i].push_back(u);
        }
    }
}

// ---------------------------------------------------------------------
// 4. Exists Solution (MNI): Use Integer IDs and Bitsets
// ---------------------------------------------------------------------
static bool exists_solution_with(const DataGraph& G, const Pattern& P,
                                 int fixVar, int fixNode,
                                 const vector<vector<int>>& domainsInit)
{
    const int k = (int)P.vlab.size();
    const int n = (int)G.vlabels.size();

    // OPTIMIZATION: Convert Labels to Integers ONCE
    vector<int> edge_ids;
    edge_ids.reserve(P.pedges.size());
    for(const auto& e : P.pedges) {
        auto it = G.label_to_id.find(e.el);
        if(it == G.label_to_id.end()) return false;
        edge_ids.push_back(it->second);
    }

    vector<int> assign(k, -1);
    vector<char> used(n, 0);
    assign[fixVar] = fixNode;
    used[fixNode] = 1;

    auto choose_var = [&]()->int{
        int best=-1, bestCnt=INT_MAX;
        for (int v=0; v<k; ++v){
            if (assign[v]!=-1) continue;
            int cnt=0;
            for (int gi : domainsInit[v]){
                if (used[gi]) continue;
                bool ok = true;
                for (size_t i = 0; i < P.pedges.size(); ++i){
                    const auto& e = P.pedges[i];
                    int lid = edge_ids[i]; // Fast ID
                    if (e.a == v && assign[e.b] != -1){
                        if (!check_edge_fast(G,e.dir,lid,gi,assign[e.b])) { ok=false; break; }
                    } else if (e.b == v && assign[e.a] != -1){
                        if (!check_edge_fast(G,e.dir,lid,assign[e.a],gi)) { ok=false; break; }
                    }
                }
                if (ok){ ++cnt; if (cnt >= bestCnt) break; }
            }
            if (cnt < bestCnt){ best=v; bestCnt=cnt; }
        }
        return best;
    };

    auto dfs = [&](auto&& self) -> bool {
        int v = choose_var();
        if (v == -1) return true;

        for (int gi : domainsInit[v]){
            if (used[gi]) continue;
            bool ok = true;
            for (size_t i = 0; i < P.pedges.size(); ++i){
                const auto& e = P.pedges[i];
                int lid = edge_ids[i];
                if (e.a == v && assign[e.b] != -1){
                    if (!check_edge_fast(G,e.dir,lid,gi,assign[e.b])) { ok=false; break; }
                } else if (e.b == v && assign[e.a] != -1){
                    if (!check_edge_fast(G,e.dir,lid,assign[e.a],gi)) { ok=false; break; }
                }
            }
            if (!ok) continue;

            assign[v]=gi; used[gi]=1;
            if (self(self)) return true;
            used[gi]=0; assign[v]=-1;
        }
        return false;
    };
    return dfs(dfs);
}

// ---------------------------------------------------------------------
// 5. Full Support Counter (SuGraMi) with Safety Cap
// ---------------------------------------------------------------------
static const long long MAX_SUPPORT_LIMIT = 2000; // Cap to prevent stalls

static long long count_total_embeddings(const DataGraph& G, const Pattern& P,
                                        const vector<vector<int>>& domains)
{
    const int k = (int)P.vlab.size();
    const int n = (int)G.vlabels.size();

    vector<int> edge_ids;
    for(const auto& e : P.pedges) {
        auto it = G.label_to_id.find(e.el);
        if(it == G.label_to_id.end()) return 0;
        edge_ids.push_back(it->second);
    }

    vector<int> assign(k, -1);
    vector<char> used(n, 0);
    long long total = 0;

    auto dfs = [&](auto&& self, int v_idx) -> void {
        if (total >= MAX_SUPPORT_LIMIT) return; // Safety Cap

        if (v_idx == k) { total++; return; }

        int v = v_idx;
        for (int gi : domains[v]){
            if (used[gi]) continue;

            bool ok = true;
            for (size_t i = 0; i < P.pedges.size(); ++i){
                const auto& e = P.pedges[i];
                if (e.a == v && e.b < v) { // neighbor e.b assigned
                     if (!check_edge_fast(G, e.dir, edge_ids[i], gi, assign[e.b])) { ok=false; break; }
                } else if (e.b == v && e.a < v) { // neighbor e.a assigned
                     if (!check_edge_fast(G, e.dir, edge_ids[i], assign[e.a], gi)) { ok=false; break; }
                }
            }
            if (!ok) continue;

            assign[v] = gi; used[gi] = 1;
            self(self, v + 1);
            used[gi] = 0; assign[v] = -1;

            if (total >= MAX_SUPPORT_LIMIT) break; // Unwind check
        }
    };
    dfs(dfs, 0);
    return total;
}

static int compute_MNI_support_exact(const DataGraph& G, const Pattern& P, int tau){
    const int k = (int)P.vlab.size();
    if (k == 0) return 0;

    vector<vector<int>> dom(k);
    for (int i=0; i<k; ++i){
        auto it = G.lab2nodes.find(P.vlab[i]);
        if (it != G.lab2nodes.end()) dom[i].assign(it->second.begin(), it->second.end());
        if ((int)dom[i].size() < tau) return 0;
    }
    filter_domains_by_local_constraints(G, P, dom);

    int support = numeric_limits<int>::max();
    for (int v=0; v<k; ++v){
        if ((int)dom[v].size() < tau) return 0;
        int count_v = 0;
        for (int u : dom[v]){
            if (exists_solution_with(G, P, v, u, dom)){
                ++count_v;
            } else {
                 // MRV Pruning: remove 'u' from dom[v] to speed up future checks
                 // For now, we rely on the check failing fast.
            }
        }
        support = min(support, count_v);
        if (support < tau) return 0;
    }
    return support;
}

// Quick seed support (2 nodes, 1 edge) from seed map
static int mni_support_seed_from_map(const DataGraph& G, const Pattern& P,
                                     const unordered_map<string,int>& seed_mni)
{
    const auto& e = P.pedges[0];
    DataGraph::EdgeTypeKey k;
    if (G.directed){
        k = { P.vlab[e.a], P.vlab[e.b], e.el, 1 };
    }else{
        const string &la = P.vlab[e.a], &lb = P.vlab[e.b];
        if (la <= lb) k = { la, lb, e.el, 0 };
        else          k = { lb, la, e.el, 0 };
    }
    auto it = seed_mni.find(encKey(k));
    return (it==seed_mni.end()? 0 : it->second);
}

static int compute_MNI_support_hybrid(const DataGraph& G, const Pattern& P, int tau,
                                      const unordered_map<string,int>& seed_mni)
{
    if (P.vlab.size()==2 && P.pedges.size()==1){
        return mni_support_seed_from_map(G, P, seed_mni);
    }
    return compute_MNI_support_exact(G, P, tau);
}

// =================== Candidate enumeration ============================

static inline bool edge_already_in_pattern(const Pattern& S,
                                           int a, int b,
                                           const string& el,
                                           int dirflag)
{
    for (const auto& e : S.pedges){
        if (e.el != el) continue;
        if (dirflag == 0 && e.dir == 0){
            if ((e.a == a && e.b == b) || (e.a == b && e.b == a)) return true;
        } else if (dirflag == 1 && e.dir == 1){
            if (e.a == a && e.b == b) return true;
        }
    }
    return false;
}

// Necessary seed lower bound: every edge-type in P must have seed MNI >= tau
static bool seed_lower_bound_ok(const DataGraph& G, const Pattern& P, int tau,
                                const unordered_map<string,int>& seed_mni)
{
    for (const auto& e : P.pedges){
        DataGraph::EdgeTypeKey k;
        if (G.directed){
            k = { P.vlab[e.a], P.vlab[e.b], e.el, 1 };
        }else{
            const string &la = P.vlab[e.a], &lb = P.vlab[e.b];
            if (la <= lb) k = { la, lb, e.el, 0 };
            else          k = { lb, la, e.el, 0 };
        }
        auto it = seed_mni.find(encKey(k));
        if (it == seed_mni.end() || it->second < tau) return false;
    }
    return true;
}

static void enumerate_candidates(const DataGraph& G,
                                 const Pattern& S,
                                 const vector<SeedInfo>& seeds,
                                 const unordered_map<string,int>& seed_mni,
                                 int tau,
                                 vector<Pattern>& out)
{
    const int k = (int)S.vlab.size();
    unordered_set<string> seen; seen.reserve(512);

    // Node-extensions
    for (int u = 0; u < k; ++u){
        for (const auto& s : seeds){
            const auto& ek = s.key;
            if (!(ek.lu == S.vlab[u] || ek.lv == S.vlab[u])) continue;

            if (ek.dirflag == 1){
                // u as source
                if (S.vlab[u] == ek.lu){
                    Pattern ext = S;
                    ext.vlab.push_back(ek.lv);
                    int newId = (int)ext.vlab.size() - 1;
                    if (!edge_already_in_pattern(S, u, newId, ek.el, 1)){
                        ext.pedges.push_back({u, newId, ek.el, 1});
                        if (seed_lower_bound_ok(G, ext, tau, seed_mni)){
                            string key = canonical_key(ext);
                            if (seen.insert(key).second) out.push_back(std::move(ext));
                        }
                    }
                }
                // u as target
                if (S.vlab[u] == ek.lv){
                    Pattern ext = S;
                    ext.vlab.push_back(ek.lu);
                    int newId = (int)ext.vlab.size() - 1;
                    if (!edge_already_in_pattern(S, newId, u, ek.el, 1)){
                        ext.pedges.push_back({newId, u, ek.el, 1});
                        if (seed_lower_bound_ok(G, ext, tau, seed_mni)){
                            string key = canonical_key(ext);
                            if (seen.insert(key).second) out.push_back(std::move(ext));
                        }
                    }
                }
            } else {
                // undirected
                if (S.vlab[u] == ek.lu){
                    Pattern ext = S;
                    ext.vlab.push_back(ek.lv);
                    int newId = (int)ext.vlab.size() - 1;
                    if (!edge_already_in_pattern(S, u, newId, ek.el, 0)){
                        ext.pedges.push_back({u, newId, ek.el, 0});
                        if (seed_lower_bound_ok(G, ext, tau, seed_mni)){
                            string key = canonical_key(ext);
                            if (seen.insert(key).second) out.push_back(std::move(ext));
                        }
                    }
                } else if (S.vlab[u] == ek.lv){
                    Pattern ext = S;
                    ext.vlab.push_back(ek.lu);
                    int newId = (int)ext.vlab.size() - 1;
                    if (!edge_already_in_pattern(S, u, newId, ek.el, 0)){
                        ext.pedges.push_back({u, newId, ek.el, 0});
                        if (seed_lower_bound_ok(G, ext, tau, seed_mni)){
                            string key = canonical_key(ext);
                            if (seen.insert(key).second) out.push_back(std::move(ext));
                        }
                    }
                }
            }
        }
    }

    // Edge-extensions
    for (int a = 0; a < k; ++a){
        for (int b = a + 1; b < k; ++b){
            for (const auto& s : seeds){
                const auto& ek = s.key;

                if (ek.dirflag == 0){
                    if (!((ek.lu==S.vlab[a] && ek.lv==S.vlab[b]) ||
                          (ek.lu==S.vlab[b] && ek.lv==S.vlab[a])))
                        continue;

                    if (edge_already_in_pattern(S, a, b, ek.el, 0)) continue;

                    Pattern ext = S;
                    ext.pedges.push_back({a, b, ek.el, 0});
                    if (seed_lower_bound_ok(G, ext, tau, seed_mni)){
                        string key = canonical_key(ext);
                        if (seen.insert(key).second) out.push_back(std::move(ext));
                    }
                } else {
                    // a->b
                    if (ek.lu==S.vlab[a] && ek.lv==S.vlab[b]){
                        if (!edge_already_in_pattern(S, a, b, ek.el, 1)){
                            Pattern ext = S;
                            ext.pedges.push_back({a, b, ek.el, 1});
                            if (seed_lower_bound_ok(G, ext, tau, seed_mni)){
                                string key = canonical_key(ext);
                                if (seen.insert(key).second) out.push_back(std::move(ext));
                            }
                        }
                    }
                    // b->a
                    if (ek.lu==S.vlab[b] && ek.lv==S.vlab[a]){
                        if (!edge_already_in_pattern(S, b, a, ek.el, 1)){
                            Pattern ext = S;
                            ext.pedges.push_back({b, a, ek.el, 1});
                            if (seed_lower_bound_ok(G, ext, tau, seed_mni)){
                                string key = canonical_key(ext);
                                if (seen.insert(key).second) out.push_back(std::move(ext));
                            }
                        }
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------
// 6. Subgraph Extension with Full Support Flag
// ---------------------------------------------------------------------
static void SUBGRAPHEXTENSION(const DataGraph& G, int tau,
                              const vector<SeedInfo>& seeds,
                              const unordered_map<string,int>& seed_mni,
                              const Pattern& S,
                              unordered_set<string>& emitted,
                              vector<Found>& out,
                              bool compute_full_support, int max_edges)
{
    if (max_edges > 0 && S.pedges.size() >= (size_t)max_edges) return;
    const string K = canonical_key(S);
    if (!emitted.insert(K).second) return;

    // 1. MNI Check (Always runs, now optimized!)
    int mni = compute_MNI_support_hybrid(G, S, tau, seed_mni);
    if (mni < tau) return;

    // 2. Full Support Check (Optional)
    long long reported_support = mni;
    if (compute_full_support) {
        const int k = (int)S.vlab.size();
        vector<vector<int>> dom(k);
        for(int i=0; i<k; ++i) {
             auto it = G.lab2nodes.find(S.vlab[i]);
             if(it != G.lab2nodes.end()) dom[i].assign(it->second.begin(), it->second.end());
        }
        filter_domains_by_local_constraints(G, S, dom);
        reported_support = count_total_embeddings(G, S, dom);

        // SuGraMi strictness:
        if (reported_support < tau) return;
    }

    out.push_back({S, reported_support});

    // 3. Extend
    vector<Pattern> cand;
    enumerate_candidates(G, S, seeds, seed_mni, tau, cand);

    for (const auto& c : cand){
        int s = compute_MNI_support_hybrid(G, c, tau, seed_mni);
        if (s >= tau){
            SUBGRAPHEXTENSION(G, tau, seeds, seed_mni, c, emitted, out, compute_full_support, max_edges);
        }
    }
}

// ======================= Prune Infrequent Elements =====================

void prune_infrequent_graph_elements(DataGraph& G, int tau) {
    if (tau <= 1) return;

    // 1. Count Frequencies
    std::unordered_map<std::string, int> node_label_counts;
    std::unordered_map<std::string, int> edge_label_counts;

    for(const auto& lab : G.vlabels) node_label_counts[lab]++;

    for(int u = 0; u < (int)G.adj.size(); ++u) {
        for(const auto& edge : G.adj[u]) {
            if (!G.directed && u > edge.first) continue;
            edge_label_counts[edge.second]++;
        }
    }

    // 2. Identify Valid Labels
    std::unordered_set<std::string> valid_node_labels;
    std::unordered_set<std::string> valid_edge_labels;
    for(const auto& kv : node_label_counts) if(kv.second >= tau) valid_node_labels.insert(kv.first);
    for(const auto& kv : edge_label_counts) if(kv.second >= tau) valid_edge_labels.insert(kv.first);

    // 3. Map Old IDs to New IDs
    std::vector<int> old_to_new(G.vlabels.size(), -1);
    std::vector<std::string> new_vlabels;
    int new_n = 0;

    for(int i = 0; i < (int)G.vlabels.size(); ++i) {
        if(valid_node_labels.count(G.vlabels[i])) {
            old_to_new[i] = new_n++;
            new_vlabels.push_back(G.vlabels[i]);
        }
    }

    if(new_n == (int)G.vlabels.size() && valid_edge_labels.size() == edge_label_counts.size()) {
        return;
    }

    // 4. Construct New Graph Structure
    std::vector<std::vector<std::pair<int, std::string>>> new_adj(new_n);
    std::vector<std::vector<std::pair<int, std::string>>> new_rev(new_n);
    std::vector<std::unordered_map<int, std::unordered_set<std::string>>> new_adj_set(new_n);
    std::vector<std::unordered_map<int, std::unordered_set<std::string>>> new_rev_set(new_n);

    for(int u = 0; u < (int)G.adj.size(); ++u) {
        if(old_to_new[u] == -1) continue;
        int new_u = old_to_new[u];

        for(const auto& edge : G.adj[u]) {
            int v = edge.first;
            const std::string& el = edge.second;
            if(old_to_new[v] != -1 && valid_edge_labels.count(el)) {
                int new_v = old_to_new[v];
                new_adj[new_u].push_back({new_v, el});
                new_adj_set[new_u][new_v].insert(el);
            }
        }
        for(const auto& edge : G.rev[u]) {
             int v = edge.first;
             const std::string& el = edge.second;
             if(old_to_new[v] != -1 && valid_edge_labels.count(el)) {
                 int new_v = old_to_new[v];
                 new_rev[new_u].push_back({new_v, el});
                 new_rev_set[new_u][new_v].insert(el);
             }
        }
    }

    // 5. Swap and Rebuild Indices
    G.vlabels = std::move(new_vlabels);
    G.adj = std::move(new_adj);
    G.rev = std::move(new_rev);
    G.adj_set = std::move(new_adj_set);
    G.rev_set = std::move(new_rev_set);

    G.build_indices();
}

// ============================ Driver ==================================

Output run_sopagrami(const DataGraph& G_in, const Params& p){
    DataGraph G = G_in;
    prune_infrequent_graph_elements(G, p.tau);

    auto seeds = compute_frequent_edge_seeds(G, p.tau);
    if (seeds.empty()) return {};

    // ASCENDING sort order for efficient pruning
    sort(seeds.begin(), seeds.end(),
         [](const SeedInfo& a, const SeedInfo& b){ return a.full < b.full; });

    auto seed_mni = build_seed_mni_map(seeds);

    int T = 1;
#ifdef _OPENMP
    T = (p.num_threads > 0 ? p.num_threads : omp_get_max_threads());
    if (T < 1) T = 1;
    omp_set_num_threads(T);
#endif

    vector<vector<Found>> locals(T);
    vector<unordered_set<string>> local_emitted(T);
    for(int t=0; t<T; ++t) {
        locals[t].reserve(4096);
        local_emitted[t].reserve(16384);
    }

#ifdef _OPENMP
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < (int)seeds.size(); ++i){
            Pattern seed;
            seed.vlab = {seeds[i].key.lu, seeds[i].key.lv};
            seed.pedges.push_back({0,1,seeds[i].key.el,seeds[i].key.dirflag});
            SUBGRAPHEXTENSION(G, p.tau, seeds, seed_mni, seed,
                              local_emitted[tid], locals[tid], p.compute_full_support,p.max_edges);
        }
    }
#else
    for (size_t i=0; i<seeds.size(); ++i){
        Pattern seed;
        seed.vlab = {seeds[i].key.lu, seeds[i].key.lv};
        seed.pedges.push_back({0,1,seeds[i].key.el,seeds[i].key.dirflag});
        SUBGRAPHEXTENSION(G, p.tau, seeds, seed_mni, seed,
                          local_emitted[0], locals[0], p.compute_full_support,p.max_edges);
    }
#endif

    Output out;
    unordered_set<string> global_emitted;
    global_emitted.reserve(65536);

    for (int t=0; t<T; ++t){
        for (auto &f : locals[t]){
            if (global_emitted.insert(canonical_key(f.pat)).second)
                out.frequent_patterns.push_back(std::move(f));
        }
    }

    stable_sort(out.frequent_patterns.begin(), out.frequent_patterns.end(),
        [](const Found& A, const Found& B){
            if (A.pat.pedges.size() != B.pat.pedges.size())
                return A.pat.pedges.size() < B.pat.pedges.size();
            return A.full_support > B.full_support;
        });

    return out;
}

} // namespace algo
