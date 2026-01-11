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

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace algo {

// =========================== DataGraph ===============================

static inline string encKey(const DataGraph::EdgeTypeKey& k){
    return k.lu + '\t' + k.lv + '\t' + k.el + '\t' + char('0' + k.dirflag);
}


void DataGraph::build_indices() {
    int n = (int)vlabels.size();

    // 1. Map String Labels to Integer IDs
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

    // 2. Initialize Vector-based Bitsets
    int num_labels = next_id;
    adj_bits.assign(num_labels, std::vector<Bitset>(n));
    rev_bits.assign(num_labels, std::vector<Bitset>(n));

    for(int l=0; l<num_labels; ++l) {
        for(int u=0; u<n; ++u) {
            adj_bits[l][u].init(n);
            rev_bits[l][u].init(n);
        }
    }

    // 3. Populate Bitsets & Lab2Nodes
    lab2nodes.clear();
    for(int i=0; i<n; ++i) lab2nodes[vlabels[i]].insert(i);

    for(int u=0; u<n; ++u) {
        for(const auto& e : adj[u]) {
            int v = e.first;
            int l_id = label_to_id[e.second]; // Safe lookup now
            
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
    build_indices(); 
}
vector<DataGraph::EdgeTypeStat> DataGraph::edge_type_counts_insertion_order() const {
    vector<EdgeTypeStat> stats;
    unordered_map<string,int> idx; idx.reserve(1<<12);
    auto idstr = [](const EdgeTypeKey& k){
        return k.lu + "\t" + k.lv + "\t" + k.el + "\t" + char('0'+k.dirflag);
    };

    int n = (int)vlabels.size();
    if(directed){
        for(int u=0; u<n; ++u){
            const string& lu = vlabels[u];
            for(auto [v,el]: adj[u]){
                const string& lv = vlabels[v];
                EdgeTypeKey k{lu, lv, el, 1};
                string s = idstr(k);
                auto it = idx.find(s);
                if(it==idx.end()){
                    idx.emplace(s, (int)stats.size());
                    stats.push_back({k,1});
                }else{
                    stats[it->second].count++;
                }
            }
        }
    }else{
        for(int u=0; u<n; ++u){
            const string& lu = vlabels[u];
            for(auto [v,el]: adj[u]){
                if(u>v) continue; // count each undirected edge once
                const string& lv = vlabels[v];
                auto a = lu<=lv? lu:lv;
                auto b = lu<=lv? lv:lu;
                EdgeTypeKey k{a,b,el,0};
                string s = idstr(k);
                auto it = idx.find(s);
                if(it==idx.end()){
                    idx.emplace(s, (int)stats.size());
                    stats.push_back({k,1});
                }else{
                    stats[it->second].count++;
                }
            }
        }
    }
    return stats;
}

// ============================ Pattern =================================

string Pattern::key() const {
    // Canonical-ish: labels as given + edges normalized (min,max,dcode,label) sorted
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

// Canonical key across permutations within equal-label groups (exact)
static string encode_with_order(const Pattern& S, const vector<int>& order){
    const int k = (int)S.vlab.size();
    vector<int> pos(k);
    for (int i=0;i<k;++i) pos[order[i]] = i;

    string s; s.reserve(k*8 + S.pedges.size()*16);
    s += "V:";
    for (int i=0;i<k;++i){ s += S.vlab[order[i]]; s += '|'; }

    vector<tuple<int,int,int,string>> es; es.reserve(S.pedges.size());
    for (const auto& e : S.pedges){
        int a = pos[e.a], b = pos[e.b];
        int dcode = 0;
        if (e.dir==1) dcode = (a<b ? 1 : 2);
        es.emplace_back(std::min(a,b), std::max(a,b), dcode, e.el);
    }
    sort(es.begin(), es.end());

    s += "E:";
    for (auto& t : es){
        s += to_string(get<0>(t))+"-"+to_string(get<1>(t))+"-"+to_string(get<2>(t))+"-"+get<3>(t)+"|";
    }
    return s;
}

static std::string canonical_key(const Pattern& S){
    const int k = (int)S.vlab.size();
    if (k<=1) return S.key();

    // group vertex indices by label (local 'groups' lives only in this function)
    std::map<std::string, std::vector<int>> groups;
    for (int i=0;i<k;++i) groups[S.vlab[i]].push_back(i);
    for (auto& kv : groups) std::sort(kv.second.begin(), kv.second.end());

    // collect labels to iterate deterministically
    std::vector<std::string> labels; labels.reserve(groups.size());
    for (auto& kv : groups) labels.push_back(kv.first);

    auto encode_with_order = [&](const std::vector<int>& order)->std::string{
        std::vector<int> pos(k);
        for (int i=0;i<k;++i) pos[order[i]] = i;

        std::string s; s.reserve(k*8 + S.pedges.size()*16);
        s += "V:";
        for (int i=0;i<k;++i){ s += S.vlab[order[i]]; s += '|'; }

        std::vector<std::tuple<int,int,int,std::string>> es; es.reserve(S.pedges.size());
        for (const auto& e : S.pedges){
            int a = pos[e.a], b = pos[e.b];
            int dcode = 0;                // 0 undirected; 1 a->b aligned to min; 2 b->a aligned
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

    // backtrack over permutations within each equal-label group
    std::function<void(int)> dfs = [&](int gi){
        if (gi == (int)labels.size()){
            std::string code = encode_with_order(current);
            if (!have || code < best){ best = std::move(code); have = true; }
            return;
        }
        const auto& g = groups[labels[gi]];
        std::vector<int> perm = g;
        do{
            size_t old = current.size();
            current.insert(current.end(), perm.begin(), perm.end());
            dfs(gi+1);
            current.resize(old);
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

// Correct undirected equal-label handling: mni = |S| (union of endpoints)
static vector<SeedInfo> compute_frequent_edge_seeds(const DataGraph& G, int tau){
    using K = DataGraph::EdgeTypeKey;

    struct AccDir { unordered_set<int> L, R; long long full=0; K key; }; // directed or undirected lu!=lv
    struct AccEq  { unordered_set<int> S;   long long full=0; K key; }; // undirected lu==lv

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
                A.L.insert(u);
                A.R.insert(v);
                A.full += 1;
            }
        }
    } else {
        for (int u=0; u<n; ++u){
            const string& lu = G.vlabels[u];
            for (auto [v, el] : G.adj[u]){
                if (u > v) continue; // each undirected edge once
                const string& lv = G.vlabels[v];

                if (lu == lv){
                    K k{lu, lv, el, 0};
                    auto &E = acc_eq[ encKey(k) ];
                    if (E.full == 0) E.key = k;
                    E.S.insert(u);
                    E.S.insert(v);
                    E.full += 1;
                } else {
                    K k;
                    int leftNode, rightNode;
                    if (lu <= lv) { k = {lu, lv, el, 0}; leftNode=u; rightNode=v; }
                    else          { k = {lv, lu, el, 0}; leftNode=v; rightNode=u; }

                    auto &A = acc_lr[ encKey(k) ];
                    if (A.full == 0) A.key = k;
                    A.L.insert(leftNode);
                    A.R.insert(rightNode);
                    A.full += 1;
                }
            }
        }
    }

    vector<SeedInfo> seeds;
    seeds.reserve(acc_lr.size() + acc_eq.size());

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

// Build O(1) seed MNI lookups
static unordered_map<string,int> build_seed_mni_map(const vector<SeedInfo>& seeds){
    unordered_map<string,int> m; m.reserve(seeds.size()*2);
    for (auto &s : seeds) m[encKey(s.key)] = s.mni;
    return m;
}

// ==================== MNI (exact, AC + MRV) ===========================


static inline bool check_edge_fast(const DataGraph& G, int dir, int label_id, int u, int v) {
    if (label_id < 0) return false; // Label doesn't exist in graph
    
    if (G.directed && dir == 1) {
        return G.adj_bits[label_id][u].test(v); // O(1) Array Access
    } else {
        return G.adj_bits[label_id][u].test(v);
    }
}

static inline bool check_edge_bitset(const DataGraph& G,
                                     const Pattern::PEdge& e,
                                     int u_graph, int v_graph) 
{
    // If pattern edge is directed u->v (dir=1)
    if (G.directed && e.dir == 1) {
        // We need to check if graph has u->v with label e.el
        // We can check: u's outgoing neighbors for v
        auto it = G.adj_el_bits[u_graph].find(e.el);
        if (it == G.adj_el_bits[u_graph].end()) return false;
        return it->second.test(v_graph);
    } 
    // If pattern edge is undirected (dir=0) OR graph is undirected
    else {
        // We checked symmetric population in load_from_lg, so adj_el_bits contains both.
        // We just check if u is connected to v via el.
        auto it = G.adj_el_bits[u_graph].find(e.el);
        if (it == G.adj_el_bits[u_graph].end()) return false;
        return it->second.test(v_graph);
    }
}
static inline bool consistent_edge_map(const DataGraph& G,
                                       const Pattern::PEdge& e,
                                       int va, int vb,
                                       int ga, int gb) {
    if (e.dir == 1) {
        if (e.a == va && e.b == vb) return G.has_edge(ga, gb, e.el);
        if (e.a == vb && e.b == va) return G.has_edge(gb, ga, e.el);
        return true; 
    } else {
        // Safe check for undirected or symmetric edges
        return G.has_edge(ga, gb, e.el) || G.has_edge(gb, ga, e.el);
    }
}
// Local AC (neighbor-existence) with scans (safe for directed + undirected)
static void filter_domains_by_local_constraints(const DataGraph& G,
                                                const Pattern& P,
                                                vector<vector<int>>& dom)
{
    const int k = (int)P.vlab.size();
    bool changed = true;
    while(changed) {
        changed = false;
        for (int v = 0; v < k; ++v){
            if (dom[v].empty()) continue;
            vector<int> keep; keep.reserve(dom[v].size());

            for (int gi : dom[v]){
                bool ok_all = true;
                for (const auto& e : P.pedges){
                    if (e.a != v && e.b != v) continue;
                    const int nb = (e.a==v? e.b : e.a);
                    const string& needLab = P.vlab[nb];
                    const string& el = e.el;

                    // Check if gi has ANY neighbor with label 'el' and node-label 'needLab'
                    // Optimization: We could check if that neighbor is also in dom[nb]
                    bool ok_edge = false;
                    const auto& neighbors = (G.directed && e.dir==1 && e.b==v) ? G.rev[gi] : G.adj[gi];
                    
                    for (auto [x, lbl] : neighbors) {
                        if (lbl == el && G.vlabels[x] == needLab){ 
                            ok_edge = true; break; 
                        }
                    }
                    if (!ok_edge){ ok_all = false; break; }
                }
                if (ok_all) keep.push_back(gi);
            }
            if (keep.size() < dom[v].size()) {
                dom[v].swap(keep);
                changed = true;
            }
        }
    }
}



// Check existence of at least ONE solution (for MNI)
static bool exists_solution_with(const DataGraph& G, const Pattern& P,
                                 int fixVar, int fixNode,
                                 const vector<vector<int>>& domainsInit)
{
    const int k = (int)P.vlab.size();
    const int n = (int)G.vlabels.size();

    // 1. OPTIMIZATION: Convert Pattern Edge Labels to Integers ONCE
    // We do this here so the recursive solver doesn't have to look up strings.
    vector<int> edge_ids;
    edge_ids.reserve(P.pedges.size());
    for(const auto& e : P.pedges) {
        auto it = G.label_to_id.find(e.el);
        if(it == G.label_to_id.end()) return false; // Label not in graph -> Impossible
        edge_ids.push_back(it->second);
    }

    vector<int> assign(k, -1);
    vector<char> used(n, 0);

    assign[fixVar] = fixNode;
    used[fixNode] = 1;

    // Helper: Select next variable (MRV heuristic)
    auto choose_var = [&]()->int{
        int best=-1, bestCnt=std::numeric_limits<int>::max();
        for (int v=0; v<k; ++v){
            if (assign[v]!=-1) continue;
            int cnt=0;
            for (int gi : domainsInit[v]){
                if (used[gi]) continue;
                bool ok = true;
                // Use index 'i' to access both the edge 'e' and its pre-calculated 'edge_ids[i]'
                for (size_t i = 0; i < P.pedges.size(); ++i){
                    const auto& e = P.pedges[i];
                    int lid = edge_ids[i]; // <--- FAST INTEGER ID

                    if (e.a == v && assign[e.b] != -1){
                        // Check: u->v (dir check inside helper)
                        if (!check_edge_fast(G, e.dir, lid, gi, assign[e.b])) { ok=false; break; }
                    } else if (e.b == v && assign[e.a] != -1){
                        // Check: v->u
                        if (!check_edge_fast(G, e.dir, lid, assign[e.a], gi)) { ok=false; break; }
                    }
                }
                if (ok){ ++cnt; if (cnt >= bestCnt) break; }
            }
            if (cnt < bestCnt){ best=v; bestCnt=cnt; }
        }
        return best;
    };


    // Using `auto&& self` style for lambda recursion is faster/cleaner in modern C++
    auto dfs = [&](auto&& self) -> bool {
        int v = choose_var();
        if (v == -1) return true; // All assigned, solution found

        for (int gi : domainsInit[v]){
            if (used[gi]) continue;

            // Verify edges with neighbors that are ALREADY assigned
            bool ok = true;
            for (size_t i = 0; i < P.pedges.size(); ++i){
                const auto& e = P.pedges[i];
                int lid = edge_ids[i]; // <--- USE PRE-FETCHED ID

                if (e.a == v && assign[e.b] != -1){
                    if (!check_edge_fast(G, e.dir, lid, gi, assign[e.b])) { ok=false; break; }
                } else if (e.b == v && assign[e.a] != -1){
                    if (!check_edge_fast(G, e.dir, lid, assign[e.a], gi)) { ok=false; break; }
                }
            }
            if (!ok) continue;

            // Assign and Recurse
            assign[v] = gi; 
            used[gi] = 1;
            
            if (self(self)) return true; // Found one! Return immediately.
            
            // Backtrack
            used[gi] = 0; 
            assign[v] = -1;
        }
        return false;
    };

    // Start recursion
    return dfs(dfs);
}
static long long count_total_embeddings(const DataGraph& G, const Pattern& P,
                                        const vector<vector<int>>& domains)
{
    const int k = (int)P.vlab.size();
    const int n = (int)G.vlabels.size();
    vector<int> assign(k, -1);
    vector<char> used(n, 0);
    long long total = 0;

    function<void(int)> dfs = [&](int v_idx){
        if (v_idx == k) {
            total++;
            return;
        }
        // Simple static ordering for full count prevents overhead of MRV re-calculation
        int v = v_idx; 

        for (int gi : domains[v]){
            if (used[gi]) continue;
            
            bool ok = true;
            for (const auto& e : P.pedges){
                // Check edges connected to already assigned nodes
                if (e.a == v && e.b < v) { // Neighbor e.b already assigned
                     if (!consistent_edge_map(G,e,e.a,e.b,gi,assign[e.b])) { ok=false; break; }
                } else if (e.b == v && e.a < v) { // Neighbor e.a already assigned
                     if (!consistent_edge_map(G,e,e.a,e.b,assign[e.a],gi)) { ok=false; break; }
                }
            }
            if (!ok) continue;

            assign[v] = gi; used[gi] = 1;
            dfs(v + 1);
            used[gi] = 0; assign[v] = -1;
        }
    };
    dfs(0);
    return total;
}
// Exact MNI: per-variable existence, with local AC
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
                // Optimization: If we only care about frequency, we could break here.
                // But for accurate MNI reporting, we keep counting.
                // if (count_v >= tau) break; 
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

// Hybrid: use O(1) seed for k=2, exact for larger patterns
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

    // Node-extensions: connect a new node to any existing node via a frequent edge-type
    for (int u = 0; u < k; ++u){
        for (const auto& s : seeds){
            const auto& ek = s.key;
            // u must match one side's label
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

    // Edge-extensions: add a new edge between existing nodes if frequent
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

// ======================= SUBGRAPHEXTENSION ============================

static void SUBGRAPHEXTENSION(const DataGraph& G, int tau,
                              const vector<SeedInfo>& seeds,
                              const unordered_map<string,int>& seed_mni,
                              const Pattern& S,
                              unordered_set<string>& emitted,
                              vector<Found>& out,
                              bool compute_full_support) // Add Flag
{
    const string K = canonical_key(S);
    if (!emitted.insert(K).second) return;

    // 1. Fast Check: Is it frequent based on MNI?
    int mni = compute_MNI_support_hybrid(G, S, tau, seed_mni);
    if (mni < tau) return;

    // 2. Exact Calculation: If requested, compute full support
    long long reported_support = mni;
    
    if (compute_full_support) {
        // Re-generate domains for the exact counter
        const int k = (int)S.vlab.size();
        vector<vector<int>> dom(k);
        for(int i=0; i<k; ++i) {
             auto it = G.lab2nodes.find(S.vlab[i]);
             if(it != G.lab2nodes.end()) dom[i].assign(it->second.begin(), it->second.end());
        }
        filter_domains_by_local_constraints(G, S, dom);
        
        reported_support = count_total_embeddings(G, S, dom);
        
        // Optional: Filter if Exact Support < Tau (SuGraMi strictness)
        // If you want to keep patterns that satisfy MNI even if exact support is low, remove this line.
        if (reported_support < tau) return; 
    }

    out.push_back({S, reported_support});

    // 3. Extend
    vector<Pattern> cand;
    enumerate_candidates(G, S, seeds, seed_mni, tau, cand);

    for (const auto& c : cand){
        // Heuristic check before recursion
        int s = compute_MNI_support_hybrid(G, c, tau, seed_mni);
        if (s >= tau){
            SUBGRAPHEXTENSION(G, tau, seeds, seed_mni, c, emitted, out, compute_full_support);
        }
    }
}
// [alg.cpp]

void prune_infrequent_graph_elements(DataGraph& G, int tau) {
    if (tau <= 1) return; // Nothing to prune

    // 1. Count Frequencies
    std::unordered_map<std::string, int> node_label_counts;
    std::unordered_map<std::string, int> edge_label_counts;
    
    for(const auto& lab : G.vlabels) node_label_counts[lab]++;
    
    // Count edge labels (careful with undirected double counting)
    for(int u = 0; u < (int)G.adj.size(); ++u) {
        for(const auto& edge : G.adj[u]) {
            // For undirected, G.adj has both u->v and v->u. 
            // We can just count everything and check threshold*2 if strictly undirected, 
            // or just count as seen. 
            // Standard approach: Count occurrences in the adjacency list.
            // If undirected, each edge appears twice, so threshold logic should align.
            // However, GraMi typically filters based on "number of edges having this label".
            if (!G.directed && u > edge.first) continue; // Count undirected once
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

    // If no nodes removed, check if edges need removal. If neither, return.
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
            
            // Keep edge only if target is valid AND edge label is valid
            if(old_to_new[v] != -1 && valid_edge_labels.count(el)) {
                int new_v = old_to_new[v];
                new_adj[new_u].push_back({new_v, el});
                new_adj_set[new_u][new_v].insert(el);
            }
        }
        
        // Do the same for rev
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

    // IMPORTANT: Call build_indices to regenerate bitsets on the smaller graph
    G.build_indices(); 
}
// ============================ Driver ==================================

Output run_sopagrami(const DataGraph& G_in, const Params& p){
    DataGraph G = G_in;
    prune_infrequent_graph_elements(G, p.tau);

    auto seeds = compute_frequent_edge_seeds(G, p.tau);
    if (seeds.empty()) return {};

    // Paper (SoGraMi) says: "edges with low frequencies will be mined first".
    // This implies ASCENDING order (Smallest -> Largest).
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
                              local_emitted[tid], locals[tid], p.compute_full_support);
        }
    }
#else
    for (size_t i=0; i<seeds.size(); ++i){
        Pattern seed;
        seed.vlab = {seeds[i].key.lu, seeds[i].key.lv};
        seed.pedges.push_back({0,1,seeds[i].key.el,seeds[i].key.dirflag});
        SUBGRAPHEXTENSION(G, p.tau, seeds, seed_mni, seed,
                          local_emitted[0], locals[0], p.compute_full_support);
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
