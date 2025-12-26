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
                    try{
                        double d = std::stod(elab);
                        long long iv = (long long)d; // toward zero
                        elab = std::to_string(iv);
                    }catch(...){}
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
            // store both ways for undirected graphs
            adj[e.v].push_back({e.u,e.el});
            adj_set[e.v][e.u].insert(e.el);
            rev[e.u].push_back({e.v,e.el});
            rev_set[e.u][e.v].insert(e.el);
        }
    }

    lab2nodes.clear();
    for(int i=0;i<n;++i) lab2nodes[vlabels[i]].insert(i);
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

static inline bool consistent_edge_map(const DataGraph& G,
                                       const Pattern::PEdge& e,
                                       int va, int vb,   // pattern endpoints
                                       int ga, int gb) { // graph nodes mapped to (va,vb)
    if (e.dir == 1) {
        if (e.a == va && e.b == vb) return G.has_edge(ga, gb, e.el);
        if (e.a == vb && e.b == va) return G.has_edge(gb, ga, e.el);
        return true;
    } else {
        if ((e.a == va && e.b == vb) || (e.a == vb && e.b == va)) {
            return G.has_edge(ga, gb, e.el) || G.has_edge(gb, ga, e.el);
        }
        return true;
    }
}

// Local AC (neighbor-existence) with scans (safe for directed + undirected)
static void filter_domains_by_local_constraints(const DataGraph& G,
                                                const Pattern& P,
                                                vector<vector<int>>& dom)
{
    const int k = (int)P.vlab.size();

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

                bool ok_edge = false;
                if (e.dir == 1){
                    if (e.a == v){
                        for (auto [x, lbl] : G.adj[gi]) {
                            if (lbl == el && G.vlabels[x] == needLab){ ok_edge = true; break; }
                        }
                    } else {
                        for (auto [x, lbl] : G.rev[gi]) {
                            if (lbl == el && G.vlabels[x] == needLab){ ok_edge = true; break; }
                        }
                    }
                } else {
                    for (auto [x, lbl] : G.adj[gi]) {
                        if (lbl == el && G.vlabels[x] == needLab){ ok_edge = true; break; }
                    }
                }

                if (!ok_edge){ ok_all = false; break; }
            }
            if (ok_all) keep.push_back(gi);
        }
        dom[v].swap(keep);
    }
}

// Existence of a full injective embedding with x_fixVar = fixNode
static bool exists_solution_with(const DataGraph& G, const Pattern& P,
                                 int fixVar, int fixNode,
                                 const vector<vector<int>>& domainsInit)
{
    const int k = (int)P.vlab.size();
    const int n = (int)G.vlabels.size();

    vector<vector<int>> dom = domainsInit; // do not mutate caller’s copy
    vector<int> assign(k, -1);
    vector<char> used(n, 0);

    assign[fixVar] = fixNode;
    used[fixNode] = 1;

    auto choose_var = [&]()->int{
        int best=-1, bestCnt=INT_MAX;
        for (int v=0; v<k; ++v){
            if (assign[v]!=-1) continue;
            int cnt=0;
            for (int gi : dom[v]){
                if (used[gi]) continue;
                bool ok = true;
                for (const auto& e : P.pedges){
                    if (e.a == v && assign[e.b] != -1){
                        ok = consistent_edge_map(G,e,e.a,e.b,gi,assign[e.b]);
                    } else if (e.b == v && assign[e.a] != -1){
                        ok = consistent_edge_map(G,e,e.a,e.b,assign[e.a],gi);
                    }
                    if (!ok) break;
                }
                if (ok){ ++cnt; if (cnt >= bestCnt) break; }
            }
            if (cnt < bestCnt){ best=v; bestCnt=cnt; }
        }
        return best;
    };

    function<bool()> dfs = [&](){
        bool done=true; for (int i=0;i<k;++i) if (assign[i]==-1){ done=false; break; }
        if (done) return true;

        int v = choose_var();
        if (v==-1) return false;

        for (int gi : dom[v]){
            if (used[gi]) continue;

            bool ok = true;
            for (const auto& e : P.pedges){
                if (e.a == v && assign[e.b] != -1){
                    ok = consistent_edge_map(G,e,e.a,e.b,gi,assign[e.b]);
                } else if (e.b == v && assign[e.a] != -1){
                    ok = consistent_edge_map(G,e,e.a,e.b,assign[e.a],gi);
                }
                if (!ok) break;
            }
            if (!ok) continue;

            assign[v]=gi; used[gi]=1;
            if (dfs()) return true;
            used[gi]=0; assign[v]=-1;
        }
        return false;
    };

    return dfs();
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
    for (int i=0;i<k;++i) if ((int)dom[i].size() < tau) return 0;

    int support = numeric_limits<int>::max();

    for (int v=0; v<k; ++v){
        vector<int> Dv = dom[v];
        int count_v = 0;

        for (int u : Dv){
            if (exists_solution_with(G, P, v, u, dom)){
                ++count_v;
            } else {
                auto &Dref = dom[v];
                auto it = find(Dref.begin(), Dref.end(), u);
                if (it != Dref.end()) Dref.erase(it);
                if ((int)Dref.size() < tau) return 0;
            }
        }

        support = min(support, count_v);
        if (support < tau) return 0;
    }

    return (support==numeric_limits<int>::max()? 0 : support);
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
                              vector<Found>& out)
{
    const string K = canonical_key(S);
    if (!emitted.insert(K).second) return;

    int mni = compute_MNI_support_hybrid(G, S, tau, seed_mni);
    if (mni < tau) return;

    out.push_back({S, (long long)mni});

    vector<Pattern> cand;
    enumerate_candidates(G, S, seeds, seed_mni, tau, cand);

    for (const auto& c : cand){
        int s = compute_MNI_support_hybrid(G, c, tau, seed_mni);
        if (s >= tau){
            SUBGRAPHEXTENSION(G, tau, seeds, seed_mni, c, emitted, out);
        }
    }
}

// ============================ Driver ==================================

Output run_sopagrami(const DataGraph& G, const Params& p){
    // 1) frequent 1-edge seeds by true MNI
    auto seeds = compute_frequent_edge_seeds(G, p.tau);
    if (seeds.empty()){
        Output out; return out;
    }

    // 2) sort seeds by full support (edge count) descending (SoGraMi ordering)
    sort(seeds.begin(), seeds.end(),
         [](const SeedInfo& a, const SeedInfo& b){ return a.full > b.full; });

    // 3) precompute O(1) seed MNI map (for fast lower bound + seed support)
    auto seed_mni = build_seed_mni_map(seeds);

    // 4) parallel expand per seed
    int T = 1;
#ifdef _OPENMP
    T = (p.num_threads > 0 ? p.num_threads : omp_get_max_threads());
    if (T < 1) T = 1;
    omp_set_num_threads(T);
#endif

    vector<vector<Found>> locals(T);
    vector<unordered_set<string>> local_emitted(T);
    for (int t=0; t<T; ++t){
        locals[t].reserve(1<<12);
        local_emitted[t].reserve(1<<14);
    }

#ifdef _OPENMP
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        auto& out_loc = locals[tid];
        auto& emit_loc = local_emitted[tid];

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < (int)seeds.size(); ++i){
            Pattern seed;
            seed.vlab = {seeds[i].key.lu, seeds[i].key.lv};
            seed.pedges.push_back({0,1,seeds[i].key.el,seeds[i].key.dirflag});
            SUBGRAPHEXTENSION(G, p.tau, seeds, seed_mni, seed, emit_loc, out_loc);
        }
    }
#else
    // Fallback single-thread
    for (size_t i=0; i<seeds.size(); ++i){
        Pattern seed;
        seed.vlab = {seeds[i].key.lu, seeds[i].key.lv};
        seed.pedges.push_back({0,1,seeds[i].key.el,seeds[i].key.dirflag});
        SUBGRAPHEXTENSION(G, p.tau, seeds, seed_mni, seed,
                          local_emitted[0], locals[0]);
    }
#endif

    // 5) merge & sort
    Output out;
    unordered_set<string> global_emitted; global_emitted.reserve(1<<16);
    for (int t=0; t<T; ++t){
        for (auto &f : locals[t]){
            string K = canonical_key(f.pat);
            if (global_emitted.insert(K).second) out.frequent_patterns.push_back(std::move(f));
        }
    }

    stable_sort(out.frequent_patterns.begin(), out.frequent_patterns.end(),
        [](const Found& A, const Found& B){
            if (A.pat.pedges.size()!=B.pat.pedges.size())
                return A.pat.pedges.size() < B.pat.pedges.size();
            return A.full_support > B.full_support;
        });

    return out;
}

} // namespace algo
