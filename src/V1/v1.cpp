#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <numeric>
#include <set>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <filesystem>
#include <random>
#include <ctime>
#include <unordered_set>
#include <metis.h>


using namespace std;
using namespace std::chrono;


// Output stream for file logging
ofstream log_file;

// Function to log to both console and file
void log_output(const string& message) {
    cout << message;
    if (log_file.is_open()) {
        log_file << message;
        log_file.flush();
    }
}

// Hash function for pair<int, int>
struct PairHash {
    size_t operator()(const pair<int, int>& p) const {
        return hash<int>{}(p.first) ^ (hash<int>{}(p.second) << 1);
    }
};

// Timer for high-precision timing
struct Timer {
    string name;
    steady_clock::time_point start;
    Timer(const string& n) : name(n), start(steady_clock::now()) {}
    ~Timer() {
        auto duration = duration_cast<microseconds>(steady_clock::now() - start).count();
        double ms = duration / 1000.0;
        ostringstream oss;
        oss << name << ": " << fixed << setprecision(7) << ms << " ms\n";
        log_output(oss.str());
    }
};

// Bipartite graph representation
struct BipartiteGraph {
    vector<vector<int>> U, V;
    vector<int> rank;
    bool preprocessed = false;

    void addEdge(size_t u, size_t v) {
        if (u >= U.size() || v >= V.size()) {
            throw runtime_error("Vertex index out of bounds: u=" + to_string(u) + ", v=" + to_string(v));
        }
        if (find(U[u].begin(), U[u].end(), v) == U[u].end()) {
            U[u].push_back(v);
            V[v].push_back(u);
        }
    }

    size_t edgeCount() const {
        size_t count = 0;
        for (const auto& neighbors : U) count += neighbors.size();
        return count;
    }

    // Validate graph indices
    void validate() const {
        for (size_t u = 0; u < U.size(); ++u) {
            for (int v : U[u]) {
                if (v < 0 || v >= static_cast<int>(V.size())) {
                    throw runtime_error("Invalid V index: v=" + to_string(v) + " for u=" + to_string(u));
                }
            }
        }
        for (size_t v = 0; v < V.size(); ++v) {
            for (int u : V[v]) {
                if (u < 0 || u >= static_cast<int>(U.size())) {
                    throw runtime_error("Invalid U index: u=" + to_string(u) + " for v=" + to_string(v));
                }
            }
        }
    }
};

// Butterfly counts structure
struct ButterflyCounts {
    long long global = 0;
    unordered_map<int, long long> per_vertex;
    unordered_map<pair<int, int>, long long, PairHash> per_edge;
    double counting_time = 0.0;

    void add_butterflies(int u1, int u2, long long cnt, int center, bool is_U_center, size_t U_size) {
        // Accumulate wedge counts; butterflies computed after aggregation
        global += cnt; // Temporary: store raw wedge count
        per_vertex[u1] += cnt;
        per_vertex[u2] += cnt;
        if (center >= 0) {
            int center_idx = is_U_center ? center : center + U_size;
            per_vertex[center_idx] += cnt;
        }
    }

    // Finalize butterfly counts by computing C(cnt, 2)
    void finalize_butterflies(size_t U_size) {
        long long new_global = 0;
        for (auto& [v, cnt] : per_vertex) {
            if (cnt > 1) {
                long long butterflies = cnt * (cnt - 1) / 2;
                new_global += butterflies;
                cnt = butterflies;
            } else {
                cnt = 0;
            }
        }
        global = new_global;
    }
};

//convert to CSR for METIS's sake
void convertBipartiteToCSR(const BipartiteGraph& graph,
                           std::vector<idx_t>& xadj,
                           std::vector<idx_t>& adjncy,
                           int& totalNodes)
{
    int uSize = graph.U.size();
    int vSize = graph.V.size();
    totalNodes = uSize + vSize;

    // Combined adjacency list for all nodes
    std::vector<std::unordered_set<idx_t>> adj(totalNodes);

    // Populate adjacency sets with bidirectional edges
    for (int u = 0; u < uSize; ++u) {
        for (int v : graph.U[u]) {
            int v_mapped = uSize + v;  // Shift v index
            adj[u].insert(v_mapped);   // u --> v
            adj[v_mapped].insert(u);   // v --> u
        }
    }

    // Build CSR format arrays
    xadj.resize(totalNodes + 1);
    adjncy.clear();

    idx_t edgeCount = 0;
    for (int i = 0; i < totalNodes; ++i) {
        xadj[i] = edgeCount;
        for (idx_t neighbor : adj[i]) {
            adjncy.push_back(neighbor);
            edgeCount++;
        }
    }
    xadj[totalNodes] = edgeCount; // Final boundary
}


// Read .txt file
BipartiteGraph readTxtToGraph(const string& filepath, unordered_map<int, string>& uLabels, unordered_map<int, string>& vLabels) {
    ifstream file(filepath);
    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + filepath);
    }

    BipartiteGraph graph;
    string line;
    getline(file, line);
    bool isWeighted = (line == "WeightedAdjacencyGraph");

    int n, m;
    file >> n >> m;
    if (n <= 0 || m < 0) {
        file.close();
        throw runtime_error("Invalid vertex or edge count: n=" + to_string(n) + ", m=" + to_string(m));
    }
    getline(file, line);

    int uSize = n / 2;
    int vSize = n - uSize;
    if (uSize <= 0 || vSize <= 0) {
        file.close();
        throw runtime_error("Invalid partition sizes: uSize=" + to_string(uSize) + ", vSize=" + to_string(vSize));
    }

    ostringstream oss;
    oss << "Reading graph: n=" << n << ", m=" << m << ", uSize=" << uSize << ", vSize=" << vSize << "\n";
    log_output(oss.str());

    try {
        graph.U.resize(uSize);
        graph.V.resize(vSize);
    } catch (const std::bad_alloc& e) {
        file.close();
        throw runtime_error("Memory allocation failed for graph: " + string(e.what()));
    }

    for (int i = 0; i < uSize; ++i) uLabels[i] = "u" + to_string(i);
    for (int i = 0; i < vSize; ++i) vLabels[i] = "v" + to_string(i);

    vector<int> edgeCounts(n);
    for (int i = 0; i < n; ++i) {
        file >> edgeCounts[i];
        if (edgeCounts[i] < 0) {
            file.close();
            throw runtime_error("Negative edge count at index " + to_string(i));
        }
    }
    getline(file, line);

    for (int u = 0; u < uSize; ++u) {
        int numEdges = (u == 0 ? edgeCounts[0] : edgeCounts[u] - edgeCounts[u - 1]);
        if (numEdges < 0 || numEdges > vSize) {
            file.close();
            throw runtime_error("Invalid number of edges for u=" + to_string(u) + ": " + to_string(numEdges));
        }
        for (int j = 0; j < numEdges; ++j) {
            int v;
            file >> v;
            if (isWeighted) {
                int weight;
                file >> weight;
            }
            oss.str(""); oss.clear();
            oss << "Edge: u=" << u << ", v=" << v << "\n";
            log_output(oss.str());
            if (v >= uSize && v < n) {
                int v_mapped = v - uSize;
                graph.addEdge(u, v_mapped);
            } else {
                log_output("Skipping invalid edge: u=" + to_string(u) + ", v=" + to_string(v) + "\n");
            }
        }
    }

    file.close();
    graph.validate();
    return graph;
}

// Read CSV file
BipartiteGraph readCSVToGraph(const string& filepath, const string& colU, const string& colV, unordered_map<int, string>& uLabels, unordered_map<int, string>& vLabels) {
    ifstream file(filepath);
    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + filepath);
    }

    BipartiteGraph graph;
    unordered_map<string, int> uMap, vMap;
    string line, token;
    int uIndex = -1, vIndex = -1;

    getline(file, line);
    stringstream headerStream(line);
    vector<string> headers;
    while (getline(headerStream, token, ',')) {
        headers.push_back(token);
    }

    for (size_t i = 0; i < headers.size(); ++i) {
        if (headers[i] == colU) uIndex = static_cast<int>(i);
        if (headers[i] == colV) vIndex = static_cast<int>(i);
    }
    if (uIndex == -1 || vIndex == -1) {
        file.close();
        throw runtime_error("Column names not found in header");
    }

    int uCounter = 0, vCounter = 0;

    while (getline(file, line)) {
        stringstream lineStream(line);
        vector<string> row;
        while (getline(lineStream, token, ',')) {
            row.push_back(token);
        }

        if (row.size() <= static_cast<size_t>(max(uIndex, vIndex))) continue;

        string uVal = row[uIndex];
        string vVal = row[vIndex];

        int uId;
        if (uMap.find(uVal) != uMap.end()) {
            uId = uMap[uVal];
        } else {
            uId = uCounter++;
            uMap[uVal] = uId;
            graph.U.emplace_back();
            uLabels[uId] = uVal;
        }

        int vId;
        if (vMap.find(vVal) != vMap.end()) {
            vId = vMap[vVal];
        } else {
            vId = vCounter++;
            vMap[vVal] = vId;
            graph.V.emplace_back();
            vLabels[vId] = vVal;
        }

        graph.addEdge(uId, vId);
    }

    file.close();
    graph.validate();
    return graph;
}

// General read graph function
BipartiteGraph readGraph(const string& filepath, const string& colU, const string& colV, unordered_map<int, string>& uLabels, unordered_map<int, string>& vLabels) {
    string extension = filesystem::path(filepath).extension().string();
    if (extension == ".csv") {
        return readCSVToGraph(filepath, colU, colV, uLabels, vLabels);
    } else if (extension == ".txt") {
        return readTxtToGraph(filepath, uLabels, vLabels);
    } else {
        throw runtime_error("Unsupported file extension: " + extension);
    }
}

// Print graph sample
void printGraphSample(const BipartiteGraph& graph, const unordered_map<int, string>& uLabels, const unordered_map<int, string>& vLabels, size_t sampleCount = 10) {
    size_t count = 0;
    for (size_t u = 0; u < graph.U.size() && count < sampleCount; ++u) {
        for (int v : graph.U[u]) {
            string uName = uLabels.count(u) ? uLabels.at(u) : to_string(u);
            string vName = vLabels.count(v) ? vLabels.at(v) : to_string(v);
            ostringstream oss;
            oss << uName << " --> " << vName << "\n";
            log_output(oss.str());
            if (++count >= sampleCount) return;
        }
    }
}

// Graph sparsification
BipartiteGraph sparsifyGraph(const BipartiteGraph& graph, double p, unordered_map<int, string>& uLabels, unordered_map<int, string>& vLabels) {
    Timer timer("Sparsification");
    if (p <= 0 || p > 1) {
        throw runtime_error("Sampling probability must be in (0, 1]");
    }

    BipartiteGraph sparseGraph;
    sparseGraph.U.resize(graph.U.size());
    sparseGraph.V.resize(graph.V.size());

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 1);

    for (size_t u = 0; u < graph.U.size(); ++u) {
        for (int v : graph.U[u]) {
            if (dis(gen) < p) {
                sparseGraph.addEdge(u, v);
            }
        }
    }

    log_output("Sparsified graph: |E| = " + to_string(sparseGraph.edgeCount()) + " (p = " + to_string(p) + ")\n");
    sparseGraph.validate();
    return sparseGraph;
}

// Preprocess graph
void preprocess_graph(BipartiteGraph& graph, const string& ranking_method = "degree") {
    Timer timer("Preprocessing");
    if (graph.preprocessed) return;

    graph.rank.resize(graph.U.size() + graph.V.size());
    iota(graph.rank.begin(), graph.rank.end(), 0);

    if (ranking_method == "degree") {
        sort(graph.rank.begin(), graph.rank.begin() + graph.U.size(), [&](int a, int b) {
            return graph.U[a].size() > graph.U[b].size();
        });
        sort(graph.rank.begin() + graph.U.size(), graph.rank.end(), [&](int a, int b) {
            return graph.V[a - graph.U.size()].size() > graph.V[b - graph.U.size()].size();
        });
    } else if (ranking_method == "side") {
        for (size_t i = 0; i < graph.U.size(); ++i) graph.rank[i] = i;
        for (size_t i = 0; i < graph.V.size(); ++i) graph.rank[graph.U.size() + i] = graph.U.size() + i;
    }

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

    // Verify sorting
    for (size_t u = 0; u < graph.U.size(); ++u) {
        for (size_t i = 1; i < graph.U[u].size(); ++i) {
            if (graph.rank[graph.U.size() + graph.U[u][i-1]] < graph.rank[graph.U.size() + graph.U[u][i]]) {
                log_output("Warning: U[" + to_string(u) + "] not sorted correctly\n");
                break;
            }
        }
    }
    for (size_t v = 0; v < graph.V.size(); ++v) {
        for (size_t i = 1; i < graph.V[v].size(); ++i) {
            if (graph.rank[graph.V[v][i-1]] < graph.rank[graph.V[v][i]]) {
                log_output("Warning: V[" + to_string(v) + "] not sorted correctly\n");
                break;
            }
        }
    }

    graph.preprocessed = true;
}

// Exact butterfly counting
ButterflyCounts exact_count(BipartiteGraph& graph, const string& aggregation = "auto") {
    Timer timer("Exact Counting");
    ButterflyCounts counts;
    auto count_start = high_resolution_clock::now();

    bool use_U = true;
    long long wedges_U = 0, wedges_V = 0;
    for (const auto& neighbors : graph.U) {
        wedges_U += (long long)neighbors.size() * (neighbors.size() - 1) / 2;
    }
    for (const auto& neighbors : graph.V) {
        wedges_V += (long long)neighbors.size() * (neighbors.size() - 1) / 2;
    }
    if (wedges_V < wedges_U) use_U = false;

    // Initialize per_vertex for all vertices
    counts.per_vertex.reserve(graph.U.size() + graph.V.size());
    for (size_t u = 0; u < graph.U.size(); ++u) {
        counts.per_vertex[static_cast<int>(u)] = 0;
    }
    for (size_t v = 0; v < graph.V.size(); ++v) {
        counts.per_vertex[static_cast<int>(graph.U.size() + v)] = 0;
    }

    unordered_map<pair<int, int>, long long, PairHash> wedge_counts;
    unordered_map<pair<int, int>, vector<pair<int, int>>, PairHash> edge_to_wedges;
    wedge_counts.reserve(graph.edgeCount());
    edge_to_wedges.reserve(graph.edgeCount());
    size_t wedge_count = 0;
    const size_t batch_size = 1000000;
    const size_t max_memory_bytes = 2ULL * 1024 * 1024 * 1024;
    const size_t wedge_size_bytes = sizeof(pair<int, int>) + sizeof(long long);
    const double batch_time = 1.0;

    auto process_batch = [&]() {
        for (auto& [key, cnt] : wedge_counts) {
            counts.add_butterflies(key.first, key.second, cnt, -1, !use_U, graph.U.size());
        }
        wedge_count = 0;
    };

    if (use_U) {
        if (graph.U.size() > 100) log_output("Processing wedge enumeration for U partition...\n");
        auto batch_start = high_resolution_clock::now();
        for (size_t u1 = 0; u1 < graph.U.size(); ++u1) {
            for (int v : graph.U[u1]) {
                for (int u2 : graph.V[v]) {
                    if (graph.rank[u2] <= graph.rank[u1]) continue;
                    wedge_counts[{u1, u2}]++;
                    edge_to_wedges[{u1, v}].emplace_back(u1, u2);
                    edge_to_wedges[{u2, v}].emplace_back(u1, u2);
                    wedge_count++;
                    if (wedge_count >= batch_size || wedge_counts.size() * wedge_size_bytes >= max_memory_bytes) {
                        process_batch();
                    }
                }
            }
            double elapsed = duration_cast<microseconds>(high_resolution_clock::now() - batch_start).count() / 1000000.0;
            if (elapsed >= batch_time || u1 == graph.U.size() - 1) {
                ostringstream oss;
                oss << "Processed " << u1 + 1 << " of " << graph.U.size() << " U vertices (" 
                    << (100.0 * (u1 + 1) / graph.U.size()) << "%)\n";
                log_output(oss.str());
                batch_start = high_resolution_clock::now();
            }
        }
    } else {
        if (graph.V.size() > 100) log_output("Processing wedge enumeration for V partition...\n");
        auto batch_start = high_resolution_clock::now();
        for (size_t v1 = 0; v1 < graph.V.size(); ++v1) {
            for (int u : graph.V[v1]) {
                for (int v2 : graph.U[u]) {
                    if (graph.rank[graph.U.size() + v2] <= graph.rank[graph.U.size() + v1]) continue;
                    wedge_counts[{v1, v2}]++;
                    edge_to_wedges[{u, v1}].emplace_back(v1, v2);
                    edge_to_wedges[{u, v2}].emplace_back(v1, v2);
                    wedge_count++;
                    if (wedge_count >= batch_size || wedge_counts.size() * wedge_size_bytes >= max_memory_bytes) {
                        process_batch();
                    }
                }
            }
            double elapsed = duration_cast<microseconds>(high_resolution_clock::now() - batch_start).count() / 1000000.0;
            if (elapsed >= batch_time || v1 == graph.V.size() - 1) {
                ostringstream oss;
                oss << "Processed " << v1 + 1 << " of " << graph.V.size() << " V vertices (" 
                    << (100.0 * (v1 + 1) / graph.V.size()) << "%)\n";
                log_output(oss.str());
                batch_start = high_resolution_clock::now();
            }
        }
    }

    // Process remaining wedges
    if (!wedge_counts.empty()) process_batch();

    // Finalize butterfly counts (convert wedge counts to butterflies)
    counts.finalize_butterflies(graph.U.size());

    // Compute per-edge butterfly counts
    if (graph.U.size() > 100) log_output("Processing per-edge butterfly counts...\n");
    counts.per_edge.reserve(graph.edgeCount());
    for (size_t u = 0; u < graph.U.size(); ++u) {
        for (int v : graph.U[u]) {
            pair<int, int> edge = {static_cast<int>(u), v};
            counts.per_edge[edge] = 0;
            if (edge_to_wedges.find(edge) != edge_to_wedges.end()) {
                for (const auto& wedge : edge_to_wedges[edge]) {
                    auto it = wedge_counts.find(wedge);
                    if (it != wedge_counts.end() && it->second >= 1) {
                        counts.per_edge[edge] += it->second * (it->second - 1) / 2;
                    }
                }
            }
        }
    }

    counts.counting_time = duration_cast<microseconds>(high_resolution_clock::now() - count_start).count() / 1000.0;
    return counts;
}
// Approximate butterfly counting
ButterflyCounts approximate_count(BipartiteGraph& graph, double sample_prob = 0.1) {
    Timer timer("Approximate Counting");
    ButterflyCounts counts;
    auto count_start = high_resolution_clock::now();

    if (sample_prob <= 0 || sample_prob > 1) {
        throw runtime_error("Sampling probability must be in (0, 1]");
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 1);

    long long total_wedges = 0;
    for (const auto& neighbors : graph.U) {
        total_wedges += (long long)neighbors.size() * (neighbors.size() - 1) / 2;
    }

    long long sample_size = max(1LL, static_cast<long long>(total_wedges * sample_prob));
    long long sampled_butterflies = 0;

    for (long long i = 0; i < sample_size; ++i) {
        size_t u1 = uniform_int_distribution<size_t>(0, graph.U.size() - 1)(gen);
        if (graph.U[u1].size() < 2) continue;
        vector<int> neighbors = graph.U[u1];
        shuffle(neighbors.begin(), neighbors.end(), gen);
        int v1 = neighbors[0];
        int v2 = neighbors[1];

        for (int u2 : graph.V[v1]) {
            if (u2 == static_cast<int>(u1)) continue;
            if (find(graph.U[u2].begin(), graph.U[u2].end(), v2) != graph.U[u2].end()) {
                sampled_butterflies++;
            }
        }
    }

    counts.global = static_cast<long long>(sampled_butterflies / (sample_prob * sample_prob));
    counts.counting_time = duration_cast<microseconds>(high_resolution_clock::now() - count_start).count() / 1000.0;

    log_output("Approximate butterfly count (p=" + to_string(sample_prob) + "): " + to_string(counts.global) + "\n");
    return counts;
}

// Tip decomposition
vector<long long> tip_decomposition(BipartiteGraph& graph) {
    Timer timer("Tip Decomposition");
    auto counts = exact_count(graph);
    priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> pq;
    set<int> peeled;
    for (const auto& [u, cnt] : counts.per_vertex) {
        if (u < 0 || u >= static_cast<int>(graph.U.size() + graph.V.size())) {
            log_output("Warning: Invalid vertex u=" + to_string(u) + " in priority queue\n");
            continue;
        }
        pq.push({cnt, u});
    }

    vector<long long> tip_numbers(graph.U.size() + graph.V.size(), 0);
    size_t iteration = 0;
    while (!pq.empty()) {
        auto [cnt, u] = pq.top(); pq.pop();
        if (peeled.find(u) != peeled.end()) continue;
        if (u < 0 || u >= static_cast<int>(graph.U.size() + graph.V.size())) {
            log_output("Warning: Invalid vertex u=" + to_string(u) + " popped from queue\n");
            continue;
        }
        peeled.insert(u);
        tip_numbers[u] = cnt;

        bool is_u_in_U = u < static_cast<int>(graph.U.size());
        vector<int> neighbors = is_u_in_U ? graph.U[u] : graph.V[u - graph.U.size()];
        for (int v : neighbors) {
            if (v < 0 || v >= static_cast<int>(graph.V.size())) {
                log_output("Warning: Invalid neighbor v=" + to_string(v) + " for u=" + to_string(u) + "\n");
                continue;
            }
            for (int u2 : graph.V[v]) {
                if (u2 < 0 || u2 >= static_cast<int>(graph.U.size())) {
                    log_output("Warning: Invalid u2=" + to_string(u2) + " for v=" + to_string(v) + "\n");
                    continue;
                }
                if (u2 == u || peeled.find(u2) != peeled.end()) continue;
                bool is_u2_in_U = u2 < static_cast<int>(graph.U.size());
                vector<int> u2_neighbors = is_u2_in_U ? graph.U[u2] : graph.V[u2 - graph.U.size()];
                vector<int> common;
                common.reserve(min(neighbors.size(), u2_neighbors.size()));
                set_intersection(
                    neighbors.begin(), neighbors.end(),
                    u2_neighbors.begin(), u2_neighbors.end(),
                    back_inserter(common)
                );

                long long delta = max(0LL, static_cast<long long>(common.size()) - 1);
                if (counts.per_vertex.find(u2) == counts.per_vertex.end()) {
                    log_output("Warning: Initializing counts.per_vertex[" + to_string(u2) + "] = 0\n");
                    counts.per_vertex[u2] = 0;
                }
                counts.per_vertex[u2] -= delta;
                pq.push({counts.per_vertex[u2], u2});

                if (++iteration % 1000 == 0) {
                    ostringstream oss;
                    oss << "Iteration " << iteration << ": u=" << u << ", v=" << v << ", u2=" << u2 
                        << ", common.size=" << common.size() << ", counts.per_vertex[" << u2 << "]=" << counts.per_vertex[u2] << "\n";
                    log_output(oss.str());
                }
            }
        }
    }
    return tip_numbers;
}

// Wing decomposition
unordered_map<pair<int, int>, long long, PairHash> wing_decomposition(BipartiteGraph& graph) {
    Timer timer("Wing Decomposition");
    auto counts = exact_count(graph);
    priority_queue<pair<long long, pair<int, int>>, vector<pair<long long, pair<int, int>>>, greater<>> pq;
    set<pair<int, int>> peeled;
    for (const auto& [e, cnt] : counts.per_edge) {
        if (e.first < 0 || e.first >= static_cast<int>(graph.U.size()) || e.second < 0 || e.second >= static_cast<int>(graph.V.size())) {
            log_output("Warning: Invalid edge (" + to_string(e.first) + "," + to_string(e.second) + ") in priority queue\n");
            continue;
        }
        pq.push({cnt, e});
    }

    unordered_map<pair<int, int>, long long, PairHash> wing_numbers;
    wing_numbers.reserve(counts.per_edge.size());
    while (!pq.empty()) {
        auto [cnt, e] = pq.top(); pq.pop();
        if (peeled.find(e) != peeled.end()) continue;
        peeled.insert(e);
        wing_numbers[e] = cnt;

        int u = e.first, v = e.second;
        if (v >= static_cast<int>(graph.V.size()) || u >= static_cast<int>(graph.U.size())) {
            log_output("Warning: Invalid edge (" + to_string(u) + "," + to_string(v) + ") in wing decomposition\n");
            continue;
        }
        for (int u2 : graph.V[v]) {
            if (u2 != u && u2 >= 0 && u2 < static_cast<int>(graph.U.size())) {
                vector<int> common;
                common.reserve(min(graph.U[u].size(), graph.U[u2].size()));
                set_intersection(
                    graph.U[u].begin(), graph.U[u].end(),
                    graph.U[u2].begin(), graph.U[u2].end(),
                    back_inserter(common)
                );
                if (common.size() >= 2) {
                    for (int v2 : common) {
                        if (v2 != v && v2 >= 0 && v2 < static_cast<int>(graph.V.size()) && peeled.find({u2, v2}) == peeled.end()) {
                            pair<int, int> edge = {u2, v2};
                            if (counts.per_edge.find(edge) == counts.per_edge.end()) {
                                log_output("Warning: Initializing counts.per_edge[(" + to_string(u2) + "," + to_string(v2) + ")] = 0\n");
                                counts.per_edge[edge] = 0;
                            }
                            counts.per_edge[edge]--;
                            pq.push({counts.per_edge[edge], edge});
                        }
                    }
                }
            }
        }
    }
    return wing_numbers;
}

// Main function
int main() {
    time_t now = time(nullptr);
    tm* ltm = localtime(&now);
    char timestamp[20];
    strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", ltm);
    string log_filename = "output_" + string(timestamp) + ".txt";
    log_file.open(log_filename);
    if (!log_file.is_open()) {
        cerr << "Failed to open log file: " << log_filename << "\n";
        return 1;
    }

    unordered_map<int, string> uLabels, vLabels;
    BipartiteGraph graph;
    string dataDir;

    if (filesystem::exists("./Data")) {
        dataDir = "./Data";
    } else if (filesystem::exists("../Data")) {
        dataDir = "../Data";
    } else {
        log_output("Data directory not found at ./Data or ../Data.\n");
        log_output("Enter path to Data directory (or press Enter to skip file input): ");
        getline(cin, dataDir);
        if (dataDir.empty() || !filesystem::exists(dataDir)) {
            dataDir = "";
            log_output("No valid Data directory provided. Only user-generated graphs are available.\n");
        }
    }

    log_output("Select input type:\n");
    if (!dataDir.empty()) {
        log_output("1. CSV file\n");
        log_output("2. TXT file\n");
    }
    log_output("Enter choice (1-" + to_string(dataDir.empty() ? 1 : 2) + "): ");
    int choice;
    cin >> choice;
    cin.ignore();

    if (!dataDir.empty() && choice == 1) {
        vector<string> csvFiles;
        for (const auto& entry : filesystem::directory_iterator(dataDir)) {
            if (entry.path().extension() == ".csv") {
                csvFiles.push_back(entry.path().filename().string());
            }
        }
        if (csvFiles.empty()) {
            log_output("No CSV files found in " + dataDir + ".\n");
            cerr << "No CSV files found in " << dataDir << ".\n";
            log_file.close();
            return 1;
        }

        log_output("\nAvailable CSV files:\n");
        for (size_t i = 0; i < csvFiles.size(); ++i) {
            log_output(to_string(i + 1) + ". " + csvFiles[i] + "\n");
        }
        log_output("Enter file number (1-" + to_string(csvFiles.size()) + "): ");
        size_t fileChoice;
        cin >> fileChoice;
        if (fileChoice < 1 || fileChoice > csvFiles.size()) {
            log_output("Invalid file choice.\n");
            cerr << "Invalid file choice.\n";
            log_file.close();
            return 1;
        }

        string filepath = dataDir + "/" + csvFiles[fileChoice - 1];
        string colU, colV;
        if (csvFiles[fileChoice - 1] == "edges.csv") {
            colU = "hero";
            colV = "comic";
        } else if (csvFiles[fileChoice - 1] == "devs.csv") {
            colU = "author_id";
            colV = "project_id";
        } else if (csvFiles[fileChoice - 1] == "moderators.csv") {
            colU = "moderator";
            colV = "subreddit";
        } else {
            log_output("Enter U column name: ");
            cin >> colU;
            log_output("Enter V column name: ");
            cin >> colV;
        }

        try {
            graph = readGraph(filepath, colU, colV, uLabels, vLabels);
        } catch (const exception& e) {
            log_output("Error reading CSV: " + string(e.what()) + "\n");
            cerr << "Error reading CSV: " << e.what() << "\n";
            log_file.close();
            return 1;
        }
    } else if (!dataDir.empty() && choice == 2) {
        vector<string> txtFiles;
        for (const auto& entry : filesystem::directory_iterator(dataDir)) {
            if (entry.path().extension() == ".txt") {
                txtFiles.push_back(entry.path().filename().string());
            }
        }
        if (txtFiles.empty()) {
            log_output("No TXT files found in " + dataDir + ".\n");
            cerr << "No TXT files found in " << dataDir << ".\n";
            log_file.close();
            return 1;
        }

        log_output("\nAvailable TXT files:\n");
        for (size_t i = 0; i < txtFiles.size(); ++i) {
            log_output(to_string(i + 1) + ". " + txtFiles[i] + "\n");
        }
        log_output("Enter file number (1-" + to_string(txtFiles.size()) + "): ");
        size_t fileChoice;
        cin >> fileChoice;
        if (fileChoice < 1 || fileChoice > txtFiles.size()) {
            log_output("Invalid file choice.\n");
            cerr << "Invalid file choice.\n";
            log_file.close();
            return 1;
        }

        string filepath = dataDir + "/" + txtFiles[fileChoice - 1];
        try {
            graph = readGraph(filepath, "", "", uLabels, vLabels);
        } catch (const exception& e) {
            log_output("Error reading TXT: " + string(e.what()) + "\n");
            cerr << "Error reading TXT: " << e.what() << "\n";
            log_file.close();
            return 1;
        }
    } else {
        log_output("Invalid input type.\n");
        cerr << "Invalid input type.\n";
        log_file.close();
        return 1;
    }

    ostringstream oss;
    oss << "Graph: |U| = " << graph.U.size() << ", |V| = " << graph.V.size() << ", |E| = " << graph.edgeCount() << "\n";
    log_output(oss.str());

    {
        Timer total_timer("Total Execution");

        log_output("\nSample edges from the graph:\n");
        printGraphSample(graph, uLabels, vLabels);

        log_output("\nUse graph sparsification? (0 for no, 1 for yes): ");
        int use_sparsification;
        cin >> use_sparsification;
        BipartiteGraph* target_graph = &graph;
        BipartiteGraph sparse_graph;
        double sparsification_p = 1.0;

        if (use_sparsification) {
            log_output("Enter sampling probability (0 to 1, e.g., 0.1): ");
            cin >> sparsification_p;
            try {
                sparse_graph = sparsifyGraph(graph, sparsification_p, uLabels, vLabels);
                target_graph = &sparse_graph;
            } catch (const exception& e) {
                log_output("Error in sparsification: " + string(e.what()) + "\n");
                cerr << "Error in sparsification: " << e.what() << "\n";
                log_file.close();
                return 1;
            }
        }

        std::vector<idx_t> xadj, adjncy;
        int totalNodes;
        convertBipartiteToCSR(*target_graph, xadj, adjncy, totalNodes);

        // Print to verify if needed
        std::cout << "CSR format for METIS:" << std::endl;
        for (int i = 0; i <= totalNodes; ++i) {
            std::cout << "xadj[" << i << "] = " << xadj[i] << std::endl;
        }
        for (size_t i = 0; i < adjncy.size(); ++i) {
            std::cout << "adjncy[" << i << "] = " << adjncy[i] << std::endl;
        }

        idx_t nvtxs = totalNodes; // number of vertices
        idx_t ncon = 1;           // number of balancing constraints
        idx_t nparts = 4;         // number of partitions you want

        std::vector<idx_t> part(nvtxs); // output partition for each node
        idx_t objval;                   // edge cut or communication volume

        int status = METIS_PartGraphKway(
            &nvtxs, &ncon, xadj.data(), adjncy.data(),
            NULL, NULL, NULL, &nparts,
            NULL, NULL, NULL,
            &objval, part.data()
        );

        if (status != METIS_OK) {
            std::cerr << "METIS partitioning failed!" << std::endl;
            return 1;
        }

        std::cout << "METIS partitioning successful. Edge cut: " << objval << std::endl;
        for (int i = 0; i < nvtxs; ++i) {
            std::cout << "Node " << i << " is in partition " << part[i] << std::endl;
        }

        preprocess_graph(*target_graph, "degree");
        auto exact_counts = exact_count(*target_graph, "auto");
        auto approx_counts = approximate_count(*target_graph, 0.1);
        auto tips = tip_decomposition(*target_graph);
        auto wings = wing_decomposition(*target_graph);

        log_output("\n===== Results =====\n");
        oss.str(""); oss.clear();
        oss << "Exact global butterfly count: " << exact_counts.global;
        if (use_sparsification) {
            oss << " (scaled: " << static_cast<long long>(exact_counts.global / (sparsification_p * sparsification_p)) << ")";
        }
        oss << "\n";
        log_output(oss.str());

        log_output("Per-vertex butterfly counts:\n");
        for (const auto& [u, cnt] : exact_counts.per_vertex) {
            if (cnt > 0) {
                oss.str(""); oss.clear();
                oss << "Vertex " << (u < static_cast<int>(graph.U.size()) ? uLabels.at(u) : vLabels.at(u - graph.U.size()))
                    << ": " << cnt;
                if (use_sparsification) {
                    oss << " (scaled: " << static_cast<long long>(cnt / (sparsification_p * sparsification_p)) << ")";
                }
                oss << "\n";
                log_output(oss.str());
            }
        }

        log_output("Per-edge butterfly counts:\n");
        size_t edge_count_limit = 100000;
        size_t edge_count_output = 0;
        for (const auto& [e, cnt] : exact_counts.per_edge) {
            if (edge_count_output >= edge_count_limit) {
                log_output("... (output truncated, too many edges)\n");
                break;
            }
            oss.str(""); oss.clear();
            oss << "Edge (" << uLabels.at(e.first) << "," << vLabels.at(e.second) << "): " << cnt;
            if (use_sparsification) {
                oss << " (scaled: " << static_cast<long long>(cnt / (sparsification_p * sparsification_p)) << ")";
            }
            oss << "\n";
            log_output(oss.str());
            edge_count_output++;
        }

        log_output("Tip numbers:\n");
        for (size_t i = 0; i < tips.size(); ++i) {
            if (tips[i] > 0) {
                oss.str(""); oss.clear();
                oss << "Vertex " << (i < graph.U.size() ? uLabels.at(i) : vLabels.at(i - graph.U.size()))
                    << ": " << tips[i];
                if (use_sparsification) {
                    oss << " (scaled: " << static_cast<long long>(tips[i] / (sparsification_p * sparsification_p)) << ")";
                }
                oss << "\n";
                log_output(oss.str());
            }
        }

        log_output("Wing numbers:\n");
        edge_count_output = 0;
        for (const auto& [e, cnt] : wings) {
            if (cnt > 0) {
                if (edge_count_output >= edge_count_limit) {
                    log_output("... (output truncated, too many edges)\n");
                    break;
                }
                oss.str(""); oss.clear();
                oss << "Edge (" << uLabels.at(e.first) << "," << vLabels.at(e.second) << "): " << cnt;
                if (use_sparsification) {
                    oss << " (scaled: " << static_cast<long long>(cnt / (sparsification_p * sparsification_p)) << ")";
                }
                oss << "\n";
                log_output(oss.str());
                edge_count_output++;
            }
        }

        oss.str(""); oss.clear();
        oss << "Exact counting time: " << fixed << setprecision(7) << exact_counts.counting_time << " ms\n";
        oss << "Approximate counting time: " << fixed << setprecision(7) << approx_counts.counting_time << " ms\n";
        log_output(oss.str());
    }

    log_file.close();
    return 0;
}
