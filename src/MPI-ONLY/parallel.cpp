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
#include <mpi.h>

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

    // Add destructor to ButterflyCounts
    ~ButterflyCounts() {
        per_vertex.clear();
        per_edge.clear();
    }

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
   
    // Read graph type
    if (!getline(file, line)) {
        file.close();
        throw runtime_error("Empty file or failed to read graph type");
    }
    bool isWeighted = (line == "WeightedAdjacencyGraph");
    log_output("Graph type: " + line + "\n");

    // Read uSize and vSize
    int uSize, vSize;
    if (!(file >> uSize)) {
        file.close();
        throw runtime_error("Failed to read uSize");
    }
    if (!(file >> vSize)) {
        file.close();
        throw runtime_error("Failed to read vSize");
    }
    if (uSize <= 0 || vSize <= 0) {
        file.close();
        throw runtime_error("Invalid vertex counts: uSize=" + to_string(uSize) + ", vSize=" + to_string(vSize));
    }
    getline(file, line); // Consume newline

    ostringstream oss;
    oss << "Reading graph: uSize=" << uSize << ", vSize=" << vSize << "\n";
    log_output(oss.str());

    // Initialize graph
    try {
        graph.U.resize(uSize);
        graph.V.resize(vSize);
    } catch (const std::bad_alloc& e) {
        file.close();
        throw runtime_error("Memory allocation failed for graph: " + string(e.what()));
    }

    // Assign labels
    for (int i = 0; i < uSize; ++i) uLabels[i] = "u" + to_string(i);
    for (int i = 0; i < vSize; ++i) vLabels[i] = "v" + to_string(i);

    // Track edges for weight assignment (if needed)
    vector<pair<int, int>> edgeList;

    // Read U vertices' neighbors (0-based indices)
    for (int u = 0; u < uSize; ++u) {
        if (!getline(file, line)) {
            file.close();
            throw runtime_error("Failed to read neighbors for u=" + to_string(u) + ": unexpected end of file");
        }
        stringstream ss(line);
        int v;
        oss.str(""); oss.clear();
        oss << "Processing U[" << u << "]: " << line << "\n";
        log_output(oss.str());
        while (ss >> v) {
            if (v < 0 || v >= vSize) {
                file.close();
                throw runtime_error("Invalid V vertex index v=" + to_string(v) + " for u=" + to_string(u));
            }
            graph.addEdge(u, v);
            edgeList.emplace_back(u, v);
        }
    }

    // Read V vertices' neighbors (0-based indices)
    for (int v = 0; v < vSize; ++v) {
        if (!getline(file, line)) {
            file.close();
            throw runtime_error("Failed to read neighbors for v=" + to_string(v) + ": unexpected end of file");
        }
        stringstream ss(line);
        int u;
        oss.str(""); oss.clear();
        oss << "Processing V[" << v << "]: " << line << "\n";
        log_output(oss.str());
        while (ss >> u) {
            if (u < 0 || u >= uSize) {
                file.close();
                throw runtime_error("Invalid U vertex index u=" + to_string(u) + " for v=" + to_string(v));
            }
            graph.addEdge(u, v); // Add edge (handles duplicates safely)
            edgeList.emplace_back(u, v);
        }
    }

    // For weighted graphs, read weights after neighbors (no strict count validation)
    if (isWeighted) {
        int edgeIndex = 0;
        while (getline(file, line)) {
            stringstream ss(line);
            int weight;
            while (ss >> weight) {
                edgeIndex++;
            }
        }
        oss.str(""); oss.clear();
        oss << "Read " << edgeIndex << " weights for " << edgeList.size() << " edges\n";
        log_output(oss.str());
        // Weights are not used, so no action needed; just log the count
    }

    // Check for extra data
    if (getline(file, line) && !line.empty()) {
        log_output("Warning: Extra data found after expected end of file: " + line + "\n");
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

void printGraphDebug(const BipartiteGraph& graph, int rank) {
    cout << "Process " << rank << " graph:" << endl;
    for (size_t u = 0; u < graph.U.size(); ++u) {
        for (int v : graph.U[u]) {
            cout << "u" << u << " --> v" << v << endl;
        }
    }
}

// Create random bipartite graph
BipartiteGraph createUserGraph(unordered_map<int, string>& uLabels, unordered_map<int, string>& vLabels) {
    BipartiteGraph graph;
    size_t uSize, vSize, numEdges;
    log_output("Enter number of U vertices: ");
    cin >> uSize;
    log_output("Enter number of V vertices: ");
    cin >> vSize;
    log_output("Enter number of edges (max " + to_string(uSize * vSize) + "): ");
    cin >> numEdges;

    if (uSize == 0 || vSize == 0) {
        throw runtime_error("Vertex counts must be positive");
    }
    if (numEdges > uSize * vSize) {
        throw runtime_error("Too many edges requested: " + to_string(numEdges));
    }

    try {
        graph.U.resize(uSize);
        graph.V.resize(vSize);
    } catch (const std::bad_alloc& e) {
        throw runtime_error("Memory allocation failed for graph: " + string(e.what()));
    }

    for (size_t i = 0; i < uSize; ++i) uLabels[i] = "u" + to_string(i);
    for (size_t i = 0; i < vSize; ++i) vLabels[i] = "v" + to_string(i);

    random_device rd;
    mt19937 gen(rd());
    vector<pair<size_t, size_t>> possibleEdges;
    possibleEdges.reserve(uSize * vSize);
    for (size_t u = 0; u < uSize; ++u) {
        for (size_t v = 0; v < vSize; ++v) {
            possibleEdges.emplace_back(u, v);
        }
    }
    shuffle(possibleEdges.begin(), possibleEdges.end(), gen);
    for (size_t i = 0; i < numEdges && i < possibleEdges.size(); ++i) {
        graph.addEdge(possibleEdges[i].first, possibleEdges[i].second);
    }

    graph.validate();
    return graph;
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

    // Initialize counts for all vertices
    for (size_t u = 0; u < graph.U.size(); ++u) {
        counts.per_vertex[static_cast<int>(u)] = 0;
    }
    for (size_t v = 0; v < graph.V.size(); ++v) {
        counts.per_vertex[static_cast<int>(graph.U.size() + v)] = 0;
    }

    // Wedge counting between U vertices
    unordered_map<pair<int, int>, long long, PairHash> wedge_counts;
    for (size_t v = 0; v < graph.V.size(); ++v) {
        const auto& neighbors = graph.V[v];
        // For each pair of U vertices connected to this V vertex
        for (size_t i = 0; i < neighbors.size(); ++i) {
            for (size_t j = i + 1; j < neighbors.size(); ++j) {
                int u1 = neighbors[i];
                int u2 = neighbors[j];
                if (u1 > u2) swap(u1, u2); // Ensure consistent ordering
                wedge_counts[{u1, u2}]++;
            }
        }
    }

    // Convert wedge counts to butterfly counts
    counts.global = 0;
    for (const auto& [pair, cnt] : wedge_counts) {
        if (cnt >= 2) {
            long long butterflies = cnt * (cnt - 1) / 2;
            counts.global += butterflies;
            // Distribute to involved vertices
            counts.per_vertex[pair.first] += butterflies;
            counts.per_vertex[pair.second] += butterflies;
        }
    }

    // Calculate per-edge counts
    counts.per_edge.clear();
    for (size_t u = 0; u < graph.U.size(); ++u) {
        for (int v : graph.U[u]) {
            long long edge_count = 0;
            // Count butterflies this edge participates in
            for (int other_u : graph.V[v]) {
                if (other_u == static_cast<int>(u)) continue;
                int u_int = static_cast<int>(u);
                auto it = wedge_counts.find({min(u_int, other_u), max(u_int, other_u)});
                if (it != wedge_counts.end() && it->second >= 2) {
                    edge_count += (it->second - 1);
                }
            }
            counts.per_edge[{static_cast<int>(u), v}] = edge_count;
        }
    }

    counts.counting_time = duration_cast<microseconds>(high_resolution_clock::now() - count_start).count() / 1000.0;
    return counts;
}

// Tip decomposition
vector<long long> tip_decomposition(BipartiteGraph& graph, BipartiteGraph& ori) {
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
            if (v < 0 || v >= static_cast<int>(ori.V.size())) {
                log_output("Warning: Invalid neighbor v=" + to_string(v) + " for u=" + to_string(u) + "\n");
                continue;
            }
            for (int u2 : graph.V[v]) {
                if (u2 < 0 || u2 >= static_cast<int>(ori.U.size())) {
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

void broadcastBipartiteGraph(BipartiteGraph& graph, int root, MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Broadcast sizes first
    size_t u_size, v_size;
    if (rank == root) {
        u_size = graph.U.size();
        v_size = graph.V.size();
    }
    MPI_Bcast(&u_size, 1, MPI_UNSIGNED_LONG, root, comm);
    MPI_Bcast(&v_size, 1, MPI_UNSIGNED_LONG, root, comm);

    // Resize vectors in non-root processes
    if (rank != root) {
        graph.U.resize(u_size);
        graph.V.resize(v_size);
    }

    // Broadcast U adjacency lists
    for (size_t u = 0; u < u_size; ++u) {
        size_t neighbor_count;
        if (rank == root) {
            neighbor_count = graph.U[u].size();
        }
        MPI_Bcast(&neighbor_count, 1, MPI_UNSIGNED_LONG, root, comm);
        
        if (rank != root) {
            graph.U[u].resize(neighbor_count);
        }
        
        if (neighbor_count > 0) {
            MPI_Bcast(graph.U[u].data(), neighbor_count, MPI_INT, root, comm);
        }
    }

    // Broadcast V adjacency lists
    for (size_t v = 0; v < v_size; ++v) {
        size_t neighbor_count;
        if (rank == root) {
            neighbor_count = graph.V[v].size();
        }
        MPI_Bcast(&neighbor_count, 1, MPI_UNSIGNED_LONG, root, comm);
        
        if (rank != root) {
            graph.V[v].resize(neighbor_count);
        }
        
        if (neighbor_count > 0) {
            MPI_Bcast(graph.V[v].data(), neighbor_count, MPI_INT, root, comm);
        }
    }
}

// Function to divide graph among processes
void divideGraph(const BipartiteGraph& full_graph, BipartiteGraph& local_graph, 
    int rank, int num_procs, unordered_map<int, string>& local_uLabels, 
    unordered_map<int, string>& local_vLabels, 
    const unordered_map<int, string>& full_uLabels, 
    const unordered_map<int, string>& full_vLabels) {

    size_t total_u = full_graph.U.size();
    size_t total_v = full_graph.V.size();

    // Calculate ranges for this process (U partition only)
    size_t u_per_proc = total_u / num_procs;
    size_t u_start = rank * u_per_proc;
    size_t u_end = (rank == num_procs - 1) ? total_u : (rank + 1) * u_per_proc;

    // Clear and resize local graph
    local_graph.U.clear();
    local_graph.V.clear();
    local_graph.U.resize(u_end - u_start);
    local_graph.V.resize(total_v); // All processes need full V for bipartite processing

    // Copy assigned U vertices
    for (size_t u = u_start; u < u_end; ++u) {
        local_graph.U[u - u_start] = full_graph.U[u];
        if (full_uLabels.count(u)) {
            local_uLabels[u - u_start] = full_uLabels.at(u);
        }
    }

    // Copy full V partition (all processes need this)
    for (size_t v = 0; v < total_v; ++v) {
        local_graph.V[v] = full_graph.V[v];
        if (full_vLabels.count(v)) {
            local_vLabels[v] = full_vLabels.at(v);
        }
    }

    // Debug print - show only this process's partition
    // cout << "Process " << rank << " partition:" << endl;
    // cout << "U vertices: ";
    // for (size_t u = 0; u < local_graph.U.size(); ++u) {
    //     cout << "u" << u + u_start << " ";
    // }
    // cout << endl << "Edges:" << endl;

    // for (size_t u = 0; u < local_graph.U.size(); ++u) {
    //     for (int v : local_graph.U[u]) {
    //         cout << "u" << u + u_start << " --> v" << v << endl;
    //     }
    // }
}

// Function to combine butterfly counts
ButterflyCounts combineButterflyCounts(const ButterflyCounts& local_counts, int root, MPI_Comm comm) {
    int rank, num_procs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_procs);
    
    ButterflyCounts global_counts;
    
    // Combine global count (butterflies are counted in all partitions they appear in)
    long long local_global = local_counts.global;
    MPI_Reduce(&local_global, &global_counts.global, 1, MPI_LONG_LONG, MPI_SUM, root, comm);
    
    // On root process: Correct the overcounting 
    if (rank == root) {
        global_counts.global /= num_procs;
    }

    // Combine per-vertex counts (butterflies are counted twice - once per U vertex)
    if (rank == root) {
        global_counts.per_vertex = local_counts.per_vertex;
        
        for (int p = 1; p < num_procs; p++) {
            // Receive size first
            size_t count;
            MPI_Recv(&count, 1, MPI_UNSIGNED_LONG, p, 0, comm, MPI_STATUS_IGNORE);
            
            // Receive pairs
            vector<pair<int, long long>> pairs(count);
            MPI_Recv(pairs.data(), count * 2, MPI_LONG_LONG, p, 0, comm, MPI_STATUS_IGNORE);
            
            // Merge into global counts (no division needed here)
            for (const auto& [v, cnt] : pairs) {
                global_counts.per_vertex[v] += cnt;
            }
        }
    } else {
        // Send per-vertex counts to root
        size_t count = local_counts.per_vertex.size();
        MPI_Send(&count, 1, MPI_UNSIGNED_LONG, root, 0, comm);
        
        vector<pair<int, long long>> pairs;
        pairs.reserve(count);
        for (const auto& [v, cnt] : local_counts.per_vertex) {
            pairs.emplace_back(v, cnt);
        }
        MPI_Send(pairs.data(), count * 2, MPI_LONG_LONG, root, 0, comm);
    }

    // Combine per-edge counts (no division needed as edges are unique to partitions)
    if (rank == root) {
        global_counts.per_edge = local_counts.per_edge; // Initialize with root's counts
        
        for (int p = 1; p < num_procs; p++) {
            // First receive the count of edges
            int edge_count;
            MPI_Recv(&edge_count, 1, MPI_INT, p, 0, comm, MPI_STATUS_IGNORE);
            
            // Then receive the edge data in a buffer
            vector<int> buffer(edge_count * 3);
            MPI_Recv(buffer.data(), edge_count * 3, MPI_INT, p, 0, comm, MPI_STATUS_IGNORE);
            
            // Process the buffer
            for (int i = 0; i < edge_count; i++) {
                int u = buffer[i*3];
                int v = buffer[i*3+1];
                long long cnt = buffer[i*3+2];
                global_counts.per_edge[{u,v}] += cnt;
            }
        }
    } else {
        // First send the count of edges
        int edge_count = local_counts.per_edge.size();
        MPI_Send(&edge_count, 1, MPI_INT, root, 0, comm);
        
        // Then pack and send the edge data
        vector<int> buffer;
        buffer.reserve(edge_count * 3);
        for (const auto& [edge, cnt] : local_counts.per_edge) {
            buffer.push_back(edge.first);
            buffer.push_back(edge.second);
            buffer.push_back(cnt);
        }
        MPI_Send(buffer.data(), edge_count * 3, MPI_INT, root, 0, comm);
    }
    return global_counts;
}

int main(int argc, char** argv) {
    // At the very beginning of main()
    high_resolution_clock::time_point program_start = high_resolution_clock::now();
    MPI_Init(&argc, &argv);
    
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    // Only rank 0 handles file I/O
    BipartiteGraph full_graph;
    unordered_map<int, string> uLabels, vLabels;
    string dataDir;
    
    if (rank == 0) {
        // Original file reading code from your main()
        time_t now = time(nullptr);
        tm* ltm = localtime(&now);
        char timestamp[20];
        strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", ltm);
        string log_filename = "output_" + string(timestamp) + ".txt";
        log_file.open(log_filename);
        
        // This reads into full_graph, uLabels, vLabels


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
            log_output((dataDir.empty() ? "1" : "3") + string(". User-generated graph\n"));
            log_output("Enter choice (1-" + to_string(dataDir.empty() ? 1 : 3) + "): ");
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
                    full_graph = readGraph(filepath, colU, colV, uLabels, vLabels);
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
                    full_graph = readGraph(filepath, "", "", uLabels, vLabels);
                } catch (const exception& e) {
                    log_output("Error reading TXT: " + string(e.what()) + "\n");
                    cerr << "Error reading TXT: " << e.what() << "\n";
                    log_file.close();
                    return 1;
                }
            } else if (dataDir.empty() ? (choice == 1) : (choice == 3)) {
                try {
                    full_graph = createUserGraph(uLabels, vLabels);
                } catch (const exception& e) {
                    log_output("Error creating user graph: " + string(e.what()) + "\n");
                    cerr << "Error creating user graph: " << e.what() << "\n";
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
            oss << "Graph: |U| = " << full_graph.U.size() << ", |V| = " << full_graph.V.size() << ", |E| = " << full_graph.edgeCount() << "\n";
            log_output(oss.str());

            {
                Timer total_timer("Total Execution");

                log_output("\nSample edges from the graph:\n");
                printGraphSample(full_graph, uLabels, vLabels);

                log_output("\nUse graph sparsification? (0 for no, 1 for yes): ");
                int use_sparsification;
                cin >> use_sparsification;
                BipartiteGraph* target_graph = &full_graph;
                BipartiteGraph sparse_graph;
                double sparsification_p = 1.0;

                if (use_sparsification) {
                    log_output("Enter sampling probability (0 to 1, e.g., 0.1): ");
                    cin >> sparsification_p;
                    try {
                        sparse_graph = sparsifyGraph(full_graph, sparsification_p, uLabels, vLabels);
                        target_graph = &sparse_graph;
                    } catch (const exception& e) {
                        log_output("Error in sparsification: " + string(e.what()) + "\n");
                        cerr << "Error in sparsification: " << e.what() << "\n";
                        log_file.close();
                        return 1;
                    }
                    full_graph = *target_graph;
                }

            }
            cout<<"Starting processing..."<<endl;
        }

    MPI_Barrier(MPI_COMM_WORLD);
    // Broadcast the graph to all processes
    broadcastBipartiteGraph(full_graph, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // Divide the graph among processes
    BipartiteGraph local_graph;
    unordered_map<int, string> local_uLabels, local_vLabels;
    divideGraph(full_graph, local_graph, rank, num_procs, 
                local_uLabels, local_vLabels, uLabels, vLabels);
  
    // Before counting, ensure no edge duplication
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "Rank " << rank << " edge count: " << local_graph.edgeCount() << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    preprocess_graph(local_graph, "degree");
    MPI_Barrier(MPI_COMM_WORLD);

    auto local_exact_counts = exact_count(local_graph, "auto");
    MPI_Barrier(MPI_COMM_WORLD);

    // After all counting is done but before output
    double local_counting_time = local_exact_counts.counting_time;
    double max_counting_time, avg_counting_time;

    // Reduce to find maximum counting time across all processes
    MPI_Reduce(&local_counting_time, &max_counting_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Reduce to find average counting time
    double sum_counting_time;
    MPI_Reduce(&local_counting_time, &sum_counting_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        avg_counting_time = sum_counting_time / num_procs;
    }

    // First ensure all processes have valid labels
    auto safe_label = [&](int id, bool is_u) -> string {
        try {
            return is_u ? local_uLabels.at(id) 
                    : local_vLabels.at(id - local_graph.U.size());
        } catch (const out_of_range&) {
            return to_string(id);  // Fallback to numeric ID
        }
    };
    
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        auto local_tips = tip_decomposition(full_graph, full_graph);
        auto local_wings = wing_decomposition(full_graph);

        // Counting Statistics (new ostringstream)
        {
            ostringstream oss;
            oss << "\n===== Counting Statistics =====" << "\n";
            oss << "Maximum process counting time: " << fixed << setprecision(7) << max_counting_time << " ms\n";
            oss << "Average process counting time: " << fixed << setprecision(7) << avg_counting_time << " ms\n";
            log_output(oss.str());
        }

        // Tip numbers (new ostringstream)
        {
            ostringstream oss;
            oss << "\n===== Root Process Tip Numbers (Sample) =====" << "\n";
            int trunc = 0;
            for (size_t i = 0; i < min(local_tips.size(), full_graph.U.size() + full_graph.V.size()); ++i) {
                if (local_tips[i] > 0 && trunc++ < 10) {
                    oss << "Vertex " << safe_label(i, i < full_graph.U.size())
                        << ": " << local_tips[i] << "\n";
                }
            }
            log_output(oss.str());
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Combine results at root
    ButterflyCounts global_exact_counts = combineButterflyCounts(local_exact_counts, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Combined Results (new ostringstream)
        ostringstream oss;
        oss << "\n===== Combined Results =====" << "\n";
        oss << "Global butterfly count: " << global_exact_counts.global << "\n";
        log_output(oss.str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Process-specific output
    for (int p = 0; p < num_procs; p++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == p) {
            ostringstream oss;
            oss << "\n===== Process " << rank << " Results =====" << "\n";
            
            // Per-vertex counts
            oss << "Per-vertex butterfly counts (sample):" << "\n";
            int vertex_count = 0;
            for (const auto& [u, cnt] : local_exact_counts.per_vertex) {
                if (cnt > 0 && vertex_count++ < 10) {
                    oss << "Vertex " << safe_label(u, u < local_graph.U.size())
                        << ": " << cnt << "\n";
                }
            }

            // Per-edge counts
            oss << "\nPer-edge butterfly counts (sample):" << "\n";
            int edge_count = 0;
            for (const auto& [e, cnt] : local_exact_counts.per_edge) {
                if (edge_count++ < 10) {
                    oss << "Edge (" << safe_label(e.first, true) << "," 
                        << safe_label(e.second, false) << "): " << cnt << "\n";
                }
            }

            oss << "\nExact counting time: " << fixed << setprecision(7) 
                << local_exact_counts.counting_time << " ms" << "\n";
            
            // Write entire process output at once
            cout << oss.str();
            if (log_file.is_open()) {
                log_file << oss.str();
                log_file.flush();
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Program timing (new ostringstream)
    if (rank == 0 && log_file.is_open()) {
        auto program_end = high_resolution_clock::now();
        auto total_duration = duration_cast<milliseconds>(program_end - program_start).count();
        ostringstream oss;
        oss << "\n===== Program Timing =====" << "\n";
        oss << "Total program execution time: " << total_duration << " ms\n";
        log_output(oss.str());
        log_file.close();
    }

    return 0;
    }