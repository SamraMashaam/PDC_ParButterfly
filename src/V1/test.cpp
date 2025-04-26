#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <string>
#include <chrono>
#include <cassert>
#include <queue>

using namespace std;

// CSR representation for bipartite graph
struct Graph {
    int n_left, n_right; // Number of vertices in left (U) and right (V) sets
    vector<int> left_adj; // Concatenated adjacency lists for left vertices
    vector<int> left_ptr; // Pointers to start of each left vertex's adjacency list
    vector<int> right_adj; // Concatenated adjacency lists for right vertices
    vector<int> right_ptr; // Pointers to start of each right vertex's adjacency list
    vector<pair<int, int>> edges; // Edge list for per-edge counting
};

// Load bipartite graph from edge list file
Graph load_graph(const string& filename) {
    Graph g;
    vector<pair<int, int>> edges;
    ifstream file(filename);
    int u, v, max_u = 0, max_v = 0;
    while (file >> u >> v) {
        edges.emplace_back(u, v);
        max_u = max(max_u, u);
        max_v = max(max_v, v);
    }
    file.close();

    g.n_left = max_u + 1;
    g.n_right = max_v + 1;

    // Count degrees for left and right vertices
    vector<int> left_deg(g.n_left, 0), right_deg(g.n_right, 0);
    for (const auto& e : edges) {
        left_deg[e.first]++;
        right_deg[e.second]++;
    }

    // Initialize CSR pointers
    g.left_ptr.resize(g.n_left + 1, 0);
    g.right_ptr.resize(g.n_right + 1, 0);
    for (int i = 0; i < g.n_left; ++i) g.left_ptr[i + 1] = g.left_ptr[i] + left_deg[i];
    for (int i = 0; i < g.n_right; ++i) g.right_ptr[i + 1] = g.right_ptr[i] + right_deg[i];

    // Fill adjacency lists
    g.left_adj.resize(g.left_ptr[g.n_left]);
    g.right_adj.resize(g.right_ptr[g.n_right]);
    vector<int> left_pos(g.n_left, 0), right_pos(g.n_right, 0);
    for (const auto& e : edges) {
        int u = e.first, v = e.second;
        g.left_adj[g.left_ptr[u] + left_pos[u]++] = v;
        g.right_adj[g.right_ptr[v] + right_pos[v]++] = u;
    }

    g.edges = edges;
    return g;
}

// Preprocessing: Approximate degree ordering (Algorithm 1)
vector<int> preprocess(const Graph& g) {
    vector<pair<int, int>> vertex_rank(g.n_left + g.n_right);
    for (int u = 0; u < g.n_left; ++u) {
        int deg = g.left_ptr[u + 1] - g.left_ptr[u];
        vertex_rank[u] = {deg, u}; // Left vertices
    }
    for (int v = 0; v < g.n_right; ++v) {
        int deg = g.right_ptr[v + 1] - g.right_ptr[v];
        vertex_rank[g.n_left + v] = {deg, v + g.n_left}; // Right vertices
    }
    // Sort by degree (approximate log-degree)
    sort(vertex_rank.begin(), vertex_rank.end());
    vector<int> rank(g.n_left + g.n_right);
    for (int i = 0; i < vertex_rank.size(); ++i) {
        rank[vertex_rank[i].second] = i;
    }
    return rank;
}

// Wedge retrieval and aggregation for per-vertex counting
vector<long long> count_per_vertex(const Graph& g, const vector<int>& rank) {
    vector<long long> vertex_butterflies(g.n_left + g.n_right, 0);
    unordered_map<pair<int, int>, int, hash<pair<int, int>>> wedge_counts;

    // Process left vertices (U)
    for (int u = 0; u < g.n_left; ++u) {
        for (int i = g.left_ptr[u]; i < g.left_ptr[u + 1]; ++i) {
            int v1 = g.left_adj[i];
            for (int j = i + 1; j < g.left_ptr[u + 1]; ++j) {
                int v2 = g.left_adj[j];
                if (rank[v1] > rank[u] && rank[v2] > rank[u]) {
                    wedge_counts[{min(v1, v2), max(v1, v2)}]++;
                }
            }
        }
    }

    // Process right vertices (V)
    for (int v = 0; v < g.n_right; ++v) {
        for (int i = g.right_ptr[v]; i < g.right_ptr[v + 1]; ++i) {
            int u1 = g.right_adj[i];
            for (int j = i + 1; j < g.right_ptr[v + 1]; ++j) {
                int u2 = g.right_adj[j];
                if (rank[u1] > rank[v + g.n_left] && rank[u2] > rank[v + g.n_left]) {
                    wedge_counts[{min(u1, u2), max(u1, u2)}]++;
                }
            }
        }
    }

    // Compute per-vertex butterfly counts
    for (const auto& [pair, count] : wedge_counts) {
        int a = pair.first, b = pair.second;
        long long butterflies = (long long)count * (count - 1) / 2;
        vertex_butterflies[a] += butterflies;
        vertex_butterflies[b] += butterflies;
    }

    return vertex_butterflies;
}

// Per-edge butterfly counting (Algorithm 4)
vector<long long> count_per_edge(const Graph& g, const vector<int>& rank) {
    vector<long long> edge_butterflies(g.edges.size(), 0);
    for (size_t e = 0; e < g.edges.size(); ++e) {
        int u = g.edges[e].first, v = g.edges[e].second;
        long long count = 0;
        // Count butterflies involving edge (u, v)
        for (int i = g.left_ptr[u]; i < g.left_ptr[u + 1]; ++i) {
            int v2 = g.left_adj[i];
            if (v2 == v) continue;
            for (int j = g.right_ptr[v]; j < g.right_ptr[v + 1]; ++j) {
                int u2 = g.right_adj[j];
                if (u2 == u) continue;
                if (rank[u2] > rank[u] && rank[v2] > rank[v + g.n_left]) {
                    count++;
                }
            }
        }
        edge_butterflies[e] = count;
    }
    return edge_butterflies;
}

// Vertex peeling for tip decomposition (Algorithm 5)
vector<int> peel_vertices(const Graph& g, vector<long long>& vertex_butterflies) {
    vector<int> peel_order;
    priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> pq;
    vector<bool> removed(g.n_left + g.n_right, false);

    // Initialize priority queue with butterfly counts
    for (int i = 0; i < g.n_left + g.n_right; ++i) {
        pq.emplace(vertex_butterflies[i], i);
    }

    while (!pq.empty()) {
        auto [butterfly_count, vertex] = pq.top();
        pq.pop();
        if (removed[vertex] || vertex_butterflies[vertex] != butterfly_count) continue;

        // Add to peel order
        peel_order.push_back(vertex);
        removed[vertex] = true;

        // Update butterfly counts for neighbors
        if (vertex < g.n_left) { // Left vertex
            for (int i = g.left_ptr[vertex]; i < g.left_ptr[vertex + 1]; ++i) {
                int v = g.left_adj[i];
                if (removed[v + g.n_left]) continue;
                vertex_butterflies[v + g.n_left] = 0; // Reset and recompute
                for (int j = g.right_ptr[v]; j < g.right_ptr[v + 1]; ++j) {
                    int u2 = g.right_adj[j];
                    if (removed[u2]) continue;
                    for (int k = j + 1; k < g.right_ptr[v + 1]; ++k) {
                        int u3 = g.right_adj[k];
                        if (removed[u3]) continue;
                        vertex_butterflies[v + g.n_left] += 1;
                    }
                }
                pq.emplace(vertex_butterflies[v + g.n_left], v + g.n_left);
            }
        } else { // Right vertex
            int v = vertex - g.n_left;
            for (int i = g.right_ptr[v]; i < g.right_ptr[v + 1]; ++i) {
                int u = g.right_adj[i];
                if (removed[u]) continue;
                vertex_butterflies[u] = 0; // Reset and recompute
                for (int j = g.left_ptr[u]; j < g.left_ptr[u + 1]; ++j) {
                    int v2 = g.left_adj[j];
                    if (removed[v2 + g.n_left]) continue;
                    for (int k = j + 1; k < g.left_ptr[u + 1]; ++k) {
                        int v3 = g.left_adj[k];
                        if (removed[v3 + g.n_left]) continue;
                        vertex_butterflies[u] += 1;
                    }
                }
                pq.emplace(vertex_butterflies[u], u);
            }
        }
    }

    return peel_order;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_graph.txt>" << endl;
        return 1;
    }

    // Load graph
    auto start = chrono::high_resolution_clock::now();
    Graph g = load_graph(argv[1]);
    auto load_time = chrono::high_resolution_clock::now();
    cout << "Graph loaded: " << g.n_left << " left vertices, " << g.n_right << " right vertices, "
         << g.edges.size() << " edges" << endl;
    cout << "Load time: " << chrono::duration<double>(load_time - start).count() << " seconds" << endl;

    // Preprocess
    start = chrono::high_resolution_clock::now();
    vector<int> rank = preprocess(g);
    auto preprocess_time = chrono::high_resolution_clock::now();
    cout << "Preprocessing time: " << chrono::duration<double>(preprocess_time - start).count() << " seconds" << endl;

    // Per-vertex butterfly counting
    start = chrono::high_resolution_clock::now();
    vector<long long> vertex_butterflies = count_per_vertex(g, rank);
    auto vertex_count_time = chrono::high_resolution_clock::now();
    cout << "Per-vertex counting time: " << chrono::duration<double>(vertex_count_time - start).count() << " seconds" << endl;

    // Per-edge butterfly counting
    start = chrono::high_resolution_clock::now();
    vector<long long> edge_butterflies = count_per_edge(g, rank);
    auto edge_count_time = chrono::high_resolution_clock::now();
    cout << "Per-edge counting time: " << chrono::duration<double>(edge_count_time - start).count() << " seconds" << endl;

    // Vertex peeling
    start = chrono::high_resolution_clock::now();
    vector<int> peel_order = peel_vertices(g, vertex_butterflies);
    auto peel_time = chrono::high_resolution_clock::now();
    cout << "Vertex peeling time: " << chrono::duration<double>(peel_time - start).count() << " seconds" << endl;

    // Output results (for validation)
    cout << "Sample per-vertex butterfly counts:" << endl;
    for (int i = 0; i < min(5, (int)vertex_butterflies.size()); ++i) {
        cout << "Vertex " << i << ": " << vertex_butterflies[i] << " butterflies" << endl;
    }
    cout << "Sample per-edge butterfly counts:" << endl;
    for (int i = 0; i < min(5, (int)edge_butterflies.size()); ++i) {
        cout << "Edge " << i << " (" << g.edges[i].first << "," << g.edges[i].second << "): "
             << edge_butterflies[i] << " butterflies" << endl;
    }
    cout << "Peel order (first 5 vertices):" << endl;
    for (int i = 0; i < min(5, (int)peel_order.size()); ++i) {
        cout << "Vertex " << peel_order[i] << endl;
    }

    return 0;
}