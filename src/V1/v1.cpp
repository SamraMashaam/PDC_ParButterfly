#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <numeric>
#include <set>
#include <chrono>
#include <iomanip>

using namespace std;
using namespace std::chrono;

// Hash function for pair<int, int> to use in unordered_map
struct PairHash {
    size_t operator()(const pair<int, int>& p) const {
        return hash<int>{}(p.first) ^ (hash<int>{}(p.second) << 1);
    }
};

// Timer to measure execution time with high precision
struct Timer {
    string name;
    steady_clock::time_point start;
    Timer(const string& n) : name(n), start(steady_clock::now()) {}
    ~Timer() {
        auto duration = duration_cast<microseconds>(steady_clock::now() - start).count();
        double ms = duration / 1000.0;
        cout << name << ": " << fixed << setprecision(7) << ms << " ms\n";
    }
};

// Bipartite graph representation
struct BipartiteGraph {
    vector<vector<int>> U, V; // Adjacency lists for partitions U and V
    vector<int> rank;         // Vertex ranks for ordering
    bool preprocessed = false;

    // Add edge (u, v) to the graph
    void addEdge(size_t u, size_t v) {
        if (u >= U.size() || v >= V.size()) {
            throw runtime_error("Vertex index out of bounds");
        }
        U[u].push_back(v);
        V[v].push_back(u);
    }
};

// Structure to store butterfly counts
struct ButterflyCounts {
    long long global = 0;
    unordered_map<int, int> per_vertex;
    unordered_map<pair<int, int>, int, PairHash> per_edge;
    double counting_time = 0.0;

    // Add butterflies for wedge endpoints and center
    // Implements Lemma 4.2 and Figure 3 for counting
    void add_butterflies(int u1, int u2, int cnt, int center, bool is_U_center, size_t U_size) {
        int butterflies = cnt * (cnt - 1) / 2;
        global += butterflies;
        per_vertex[u1] += butterflies;
        per_vertex[u2] += butterflies;
        if (cnt > 1) {
            int center_idx = is_U_center ? center + U_size : center;
            per_vertex[center_idx] += butterflies; // Center gets butterfly count
        }
    }
};

// Preprocess graph: assign ranks and sort neighbors
// Implements Algorithm 1 (Section 4.1)
void preprocess_graph(BipartiteGraph& graph, const string& ranking_method = "degree") {
    Timer timer("Preprocessing");
    if (graph.preprocessed) return;

    graph.rank.resize(graph.U.size() + graph.V.size());
    iota(graph.rank.begin(), graph.rank.end(), 0);

    if (ranking_method == "degree") {
        // Sort U vertices by decreasing degree
        sort(graph.rank.begin(), graph.rank.begin() + graph.U.size(), [&](int a, int b) {
            return graph.U[a].size() > graph.U[b].size();
        });
        // Sort V vertices by decreasing degree
        sort(graph.rank.begin() + graph.U.size(), graph.rank.end(), [&](int a, int b) {
            return graph.V[a - graph.U.size()].size() > graph.V[b - graph.U.size()].size();
        });
    } else if (ranking_method == "side") {
        // Side order: U vertices 0 to |U|-1, V vertices |U| to |U|+|V|-1
        for (size_t i = 0; i < graph.U.size(); ++i) graph.rank[i] = i;
        for (size_t i = 0; i < graph.V.size(); ++i) graph.rank[graph.U.size() + i] = graph.U.size() + i;
    }

    // Sort neighbors by decreasing rank
    for (auto& neighbors : graph.U) {
        sort(neighbors.begin(), neighbors.end(), [&](int a, int b) {
            return graph.rank[graph.U.size() + a] > graph.rank[graph.U.size() + b];
        });
    }
    for (auto& neighbors : graph.V) {
        sort(neighbors.begin(), neighbors.end(), [&](int a, int b) {
            return graph.rank[a] > graph.rank[b];
        });
    }
    graph.preprocessed = true;
}

// Exact butterfly counting
// Implements Figure 2 (Section 3.1)
ButterflyCounts exact_count(BipartiteGraph& graph, const string& aggregation = "auto") {
    Timer timer("Counting");
    ButterflyCounts counts;
    vector<tuple<int, int, int>> wedges; // (endpoint1, endpoint2, center)

    // Choose bipartition with fewer wedges (Section 4.3.1)
    bool use_U = true;
    long long wedges_U = 0, wedges_V = 0;
    for (const auto& neighbors : graph.U) wedges_U += (long long)neighbors.size() * (neighbors.size() - 1) / 2;
    for (const auto& neighbors : graph.V) wedges_V += (long long)neighbors.size() * (neighbors.size() - 1) / 2;
    if (wedges_V < wedges_U) use_U = false;

    auto count_start = high_resolution_clock::now();

    if (use_U) {
        // Enumerate wedges (u1, v, u2) where rank[u2] > rank[u1]
        for (size_t u1 = 0; u1 < graph.U.size(); ++u1) {
            for (int v : graph.U[u1]) {
                for (int u2 : graph.V[v]) {
                    if (graph.rank[u2] > graph.rank[u1]) {
                        wedges.emplace_back(u1, u2, v);
                    }
                }
            }
        }
    } else {
        // Enumerate wedges (v1, u, v2) where rank[v2] > rank[v1]
        for (size_t v1 = 0; v1 < graph.V.size(); ++v1) {
            for (int u : graph.V[v1]) {
                for (int v2 : graph.U[u]) {
                    if (graph.rank[graph.U.size() + v2] > graph.rank[graph.U.size() + v1]) {
                        wedges.emplace_back(v1, v2, u);
                    }
                }
            }
        }
    }

    // Aggregate wedges
    if (wedges.empty()) {
        counts.counting_time = duration_cast<microseconds>(high_resolution_clock::now() - count_start).count() / 1000.0;
        return counts;
    }

    if (aggregation == "hash" || (aggregation == "auto" && wedges.size() < 100000)) {
        unordered_map<pair<int, int>, int, PairHash> wedge_counts;
        for (const auto& [e1, e2, center] : wedges) {
            wedge_counts[{e1, e2}]++;
        }
        for (const auto& [key, cnt] : wedge_counts) {
            counts.add_butterflies(key.first, key.second, cnt, get<2>(wedges[0]), !use_U, graph.U.size());
        }
    } else if (aggregation == "histogram" && wedges.size() < 1000000) {
        // Histogram aggregation for medium-sized graphs
        vector<int> wedge_counts(wedges.size());
        unordered_map<pair<int, int>, int, PairHash> wedge_indices;
        size_t idx = 0;
        for (const auto& [e1, e2, center] : wedges) {
            auto key = pair{e1, e2};
            if (wedge_indices.find(key) == wedge_indices.end()) {
                wedge_indices[key] = idx++;
            }
            wedge_counts[wedge_indices[key]]++;
        }
        for (const auto& [key, i] : wedge_indices) {
            counts.add_butterflies(key.first, key.second, wedge_counts[i], get<2>(wedges[0]), !use_U, graph.U.size());
        }
    } else {
        // Sort-based aggregation
        sort(wedges.begin(), wedges.end(), [](const auto& a, const auto& b) {
            return tie(get<0>(a), get<1>(a)) < tie(get<0>(b), get<1>(b));
        });
        for (size_t i = 0; i < wedges.size(); ) {
            size_t j = i;
            while (j < wedges.size() && get<0>(wedges[j]) == get<0>(wedges[i]) &&
                   get<1>(wedges[j]) == get<1>(wedges[i])) j++;
            counts.add_butterflies(get<0>(wedges[i]), get<1>(wedges[i]), j - i, get<2>(wedges[i]), !use_U, graph.U.size());
            i = j;
        }
    }

    // Compute per-edge counts (Lemma 4.2)
    for (size_t u = 0; u < graph.U.size(); ++u) {
        for (int v : graph.U[u]) {
            counts.per_edge[{static_cast<int>(u), v}] = 0; // Initialize
            for (size_t u2 = 0; u2 < graph.U.size(); ++u2) {
                if (u2 == u) continue;
                vector<int> common;
                set_intersection(
                    graph.U[u].begin(), graph.U[u].end(),
                    graph.U[u2].begin(), graph.U[u2].end(),
                    back_inserter(common)
                );
                if (common.size() >= 2) { // At least two common V vertices form a butterfly
                    for (int v2 : common) {
                        if (v2 == v) {
                            counts.per_edge[{static_cast<int>(u), v}]++;
                        }
                    }
                }
            }
        }
    }

    counts.counting_time = duration_cast<microseconds>(high_resolution_clock::now() - count_start).count() / 1000.0;
    return counts;
}

// Tip decomposition (vertex peeling)
// Implements Figure 4 and Algorithm 5 (Section 3.2)
vector<int> tip_decomposition(BipartiteGraph& graph) {
    Timer timer("Vertex Peeling");
    auto counts = exact_count(graph);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;
    set<int> peeled;
    for (const auto& [u, cnt] : counts.per_vertex) {
        pq.push({cnt, u});
    }

    vector<int> tip_numbers(graph.U.size() + graph.V.size(), 0);
    while (!pq.empty()) {
        auto [cnt, u] = pq.top(); pq.pop();
        if (peeled.find(u) != peeled.end()) continue;
        peeled.insert(u);
        tip_numbers[u] = cnt;

        // Update counts for unpeeled vertices
        vector<int> neighbors = u < static_cast<int>(graph.U.size()) ? graph.U[u] : graph.V[u - graph.U.size()];
        for (int v : neighbors) {
            for (int u2 : graph.V[v]) {
                if (u2 != u && peeled.find(u2) == peeled.end()) {
                    vector<int> common;
                    set_intersection(
                        (u2 < static_cast<int>(graph.U.size()) ? graph.U[u2] : graph.V[u2 - graph.U.size()]).begin(),
                        (u2 < static_cast<int>(graph.U.size()) ? graph.U[u2] : graph.V[u2 - graph.U.size()]).end(),
                        neighbors.begin(), neighbors.end(),
                        back_inserter(common)
                    );
                    counts.per_vertex[u2] -= max(0, static_cast<int>(common.size()) - 1);
                    pq.push({counts.per_vertex[u2], u2});
                }
            }
        }
    }
    return tip_numbers;
}

// Wing decomposition (edge peeling)
// Implements Figure 4 and Algorithm 6 (Section 3.2)
unordered_map<pair<int, int>, int, PairHash> wing_decomposition(BipartiteGraph& graph) {
    Timer timer("Edge Peeling");
    auto counts = exact_count(graph);
    priority_queue<pair<int, pair<int, int>>, vector<pair<int, pair<int, int>>>, greater<>> pq;
    set<pair<int, int>> peeled;
    for (const auto& [e, cnt] : counts.per_edge) {
        pq.push({cnt, e});
    }

    unordered_map<pair<int, int>, int, PairHash> wing_numbers;
    while (!pq.empty()) {
        auto [cnt, e] = pq.top(); pq.pop();
        if (peeled.find(e) != peeled.end()) continue;
        peeled.insert(e);
        wing_numbers[e] = cnt;

        // Update counts for unpeeled edges
        int u = e.first, v = e.second;
        for (int u2 : graph.V[v]) {
            if (u2 != u) {
                vector<int> common;
                set_intersection(
                    graph.U[u].begin(), graph.U[u].end(),
                    graph.U[u2].begin(), graph.U[u2].end(),
                    back_inserter(common)
                );
                if (common.size() >= 2) {
                    for (int v2 : common) {
                        if (v2 != v && peeled.find({u2, v2}) == peeled.end()) {
                            counts.per_edge[{u2, v2}]--;
                            pq.push({counts.per_edge[{u2, v2}], {u2, v2}});
                        }
                    }
                }
            }
        }
    }
    return wing_numbers;
}

// Main function with sample graph
int main() {
    BipartiteGraph graph;
    graph.U.resize(3); // Vertices u0, u1, u2
    graph.V.resize(3); // Vertices v0, v1, v2

    // Add edges for a small bipartite graph
    graph.addEdge(0, 0); // u0-v0
    graph.addEdge(0, 1); // u0-v1
    graph.addEdge(1, 0); // u1-v0
    graph.addEdge(1, 1); // u1-v1
    graph.addEdge(2, 1); // u2-v1
    graph.addEdge(2, 2); // u2-v2

    {
        Timer total_timer("Total Execution");

        // Preprocess with degree ordering
        preprocess_graph(graph, "degree");

        // Compute butterfly counts
        auto counts = exact_count(graph, "hash");

        // Perform tip and wing decompositions
        auto tips = tip_decomposition(graph);
        auto wings = wing_decomposition(graph);

        // Print results
        cout << "\n===== Results =====\n";
        cout << "Global butterfly count: " << counts.global << "\n";
        cout << "Per-vertex butterfly counts:\n";
        for (const auto& [u, cnt] : counts.per_vertex) {
            if (cnt > 0) {
                cout << "Vertex " << (u < static_cast<int>(graph.U.size()) ? "u" : "v")
                     << (u < static_cast<int>(graph.U.size()) ? u : u - graph.U.size())
                     << ": " << cnt << "\n";
            }
        }
        cout << "Per-edge butterfly counts:\n";
        for (const auto& [e, cnt] : counts.per_edge) {
            cout << "Edge (u" << e.first << ",v" << e.second << "): " << cnt << "\n";
        }
        cout << "Tip numbers:\n";
        for (size_t i = 0; i < tips.size(); ++i) {
            if (tips[i] > 0) {
                cout << "Vertex " << (i < graph.U.size() ? "u" : "v")
                     << (i < graph.U.size() ? i : i - graph.U.size())
                     << ": " << tips[i] << "\n";
            }
        }
        cout << "Wing numbers:\n";
        for (const auto& [e, cnt] : wings) {
            if (cnt > 0) {
                cout << "Edge (u" << e.first << ",v" << e.second << "): " << cnt << "\n";
            }
        }
        cout << "Counting time: " << fixed << setprecision(7) << counts.counting_time << " ms\n";
    }

    return 0;
}