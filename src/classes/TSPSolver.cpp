#include "../headers/TSPSolver.h"


void TSPSolver::printClustersHelper() {
    cout << endl << "Greater numbers of clusters and iterations will generally lead to a less costly tour." << endl;
    cout << "Running times may vary." << endl << endl;
    this_thread::sleep_for(chrono::milliseconds(3000));
    cout << "Here are some suggestions for those numbers:" << endl;
    cout << "Clusters: 5, 10, 20, 50, 100, 200, ..." << endl;
    cout << "Iterations: 50, 100, 200, 500, 1000, 2000, ..." << endl << endl;
    this_thread::sleep_for(chrono::milliseconds(2000));
    cout << "However, the choices for these values are completely up you." << endl << endl;
    this_thread::sleep_for(chrono::milliseconds(2000));
    cout << "You are encouraged to try multiple combinations of inputs." << endl << endl;
    this_thread::sleep_for(chrono::milliseconds(2000));
    cout << "Larger numbers will typically take more time." << endl << endl;
    this_thread::sleep_for(chrono::milliseconds(2000));
}


/**
 * @brief Calculates the Haversine distance between two vertices.
 * @param v1 Pointer to the first vertex.
 * @param v2 Pointer to the second vertex.
 * @return The Haversine distance in kilometers.
 * @details This function uses the Haversine formula to calculate the great-circle distance between two points on the Earth.
 *   - Time Complexity: O(1)
 *   - Space Complexity: O(1)
 */

double TSPSolver::haversineDistance(Vertex<int> *v1, Vertex<int> *v2) {
    constexpr double R = 6371.0; // Earth radius in kilometers

    double lat1 = v1->getLatitude();
    double lon1 = v1->getLongitude();
    double lat2 = v2->getLatitude();
    double lon2 = v2->getLongitude();

    lat1 *= M_PI / 180.0;
    lon1 *= M_PI / 180.0;
    lat2 *= M_PI / 180.0;
    lon2 *= M_PI / 180.0;

    double dlat = lat2 - lat1;
    double dlon = lon2 - lon1;

    double a = sin(dlat / 2) * sin(dlat / 2) +
               cos(lat1) * cos(lat2) *
               sin(dlon / 2) * sin(dlon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    double distance = R * c;

    return distance;
}

/**
 * @brief Solves the TSP using a backtracking approach.
 * @param g Pointer to the graph.
 * @param current Pointer to the current vertex.
 * @param current_cost The current cost of the path.
 * @param num_visited The number of visited vertices.
 * @param min_cost Reference to the minimum cost found.
 * @param tsp_path Reference to the vector storing the TSP path.
 * @param current_path Reference to the vector storing the current path.
 * @details This function recursively explores all possible paths to solve the TSP.
 *   - Time Complexity: O(n!), where n is the number of vertices.
 *   - Space Complexity: O(n), for the recursion stack and path storage.
 */

void TSPSolver::tspBacktrack(Graph<int>* g, Vertex<int>* current, double current_cost, int num_visited, double &min_cost, vector<Vertex<int> *> &tsp_path, vector<Vertex<int> *> &current_path) {
    current_path.push_back(current);

    if (num_visited == g->getNumVertex()) {
        for (Edge<int> *e : current->getAdj()) {
            if (e->getDest()->getInfo() == 0) { // Check if there is an edge back to the start
                double cost = current_cost + e->getWeight();
                if (cost < min_cost) {
                    min_cost = cost;
                    tsp_path = current_path; // Update tsp_path with the current best path
                    tsp_path.push_back(g->findVertex(0)); // Complete the cycle
                }
                current_path.pop_back();
                return;
            }
        }
    }

    for (Edge<int> *e : current->getAdj()) {
        Vertex<int> *w = e->getDest();
        if (!w->isVisited() && e->getWeight() != 0) { // Ensure the edge has a valid weight
            w->setVisited(true);
            tspBacktrack(g, w, current_cost + e->getWeight(), num_visited + 1, min_cost, tsp_path, current_path);
            w->setVisited(false);
        }
    }

    current_path.pop_back();
}

/**
 * @brief Solves the TSP using a brute force approach.
 * @param g Pointer to the graph.
 * @param tsp_path Reference to the vector storing the TSP path.
 * @return The minimum cost of the TSP path.
 * @details This function initializes the graph and calls the backtracking function to solve the TSP.
 *   - Time Complexity: O(n!), where n is the number of vertices.
 *   - Space Complexity: O(n), for path storage.
 */

double TSPSolver::tspBruteForce(Graph<int>* g, vector<Vertex<int> *> &tsp_path) {
    for (auto v: g->getVertexSet()) {
        v->setVisited(false);
        v->setPath(nullptr);
    }

    double min_cost = INF;
    vector<Vertex<int> *> current_path;
    Vertex<int> *source = g->findVertex(0);
    source->setVisited(true);
    tspBacktrack(g, source, 0, 1, min_cost, tsp_path, current_path);

    return min_cost;
}

/**
 * @brief Calculates the TSP using brute force and prints the result.
 * @param g Pointer to the graph.
 * @details This function calculates the TSP path and cost using the brute force method, and prints the results including elapsed time.
 *   - Time Complexity: O(n!), where n is the number of vertices.
 *   - Space Complexity: O(n), for path storage.
 */

void TSPSolver::calculateTSP(Graph<int>* g) {
    vector<Vertex<int> *> tsp_path;

    auto start = std::chrono::high_resolution_clock::now();
    double cost = tspBruteForce(g, tsp_path);
    auto end = std::chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    cout << "TSP Path : ";
    for (size_t i = 0; i < tsp_path.size(); i++) {
        cout << tsp_path[i]->getInfo() << (i == tsp_path.size() - 1 ? "\n" : " -> ");
    }
    cout << '\n';
    cout << "Cost: " << cost << '\n';
    cout << "Elapsed Time: " << duration.count() << " s\n\n";
}

/**
 * @brief Performs a preorder traversal of the MST to generate the TSP path.
 * @param g Pointer to the graph.
 * @param current Pointer to the current vertex.
 * @param result Reference to the vector storing the TSP path.
 * @param cost Reference to the total cost of the path.
 * @param prev Reference to the previous vertex in the path.
 * @details This function recursively traverses the MST in preorder to generate a TSP path.
 *   - Time Complexity: O(V + E), where V is the number of vertices and E is the number of edges.
 *   - Space Complexity: O(V), for path storage and recursion stack.
 */

void TSPSolver::preorderMST(Graph<int>* g, Vertex<int>* current, std::vector<Vertex<int>*> &result, double &cost, Vertex<int>* &prev){
    Vertex<int>* v = g->findVertex(current->getInfo());
    result.push_back(v);
    current->setVisited(true);
    bool flag = true;
    for (Edge<int>* e: current->getAdj()) {
        if (flag && !e->getDest()->isVisited()) {
            cost += e->getWeight();
            prev = e->getDest();
        }
        else if (!e->getDest()->isVisited()) {
            double tmp = 0;
            for (Edge<int> * edg :prev->getAdj() )
                if (edg->getDest()->getInfo() == e->getDest()->getInfo()){
                    tmp = edg->getWeight();
                }

            if (tmp == 0) {
                Vertex<int>* dest = e->getDest();
                if (prev != nullptr && dest != nullptr) {
                    double tmp = haversineDistance(prev, dest);
                    cost += tmp;
                }
            } else {
                cost += tmp;
            }
            prev = e->getDest();
        }
        flag = false;
        Vertex<int>* w = e->getDest();
        if (!w->isVisited()) {
            preorderMST(g, w, result, cost, prev);
        }
    }
}

/**
 * @brief Implements Prim's algorithm to generate an MST and uses preorder traversal for TSP.
 * @param g Pointer to the graph.
 * @param source Pointer to the source vertex.
 * @param result Reference to the vector storing the TSP path.
 * @param cost Reference to the total cost of the path.
 * @details This function uses Prim's algorithm to generate a Minimum Spanning Tree (MST) and then performs a preorder traversal to estimate a TSP path.
 *   - Time Complexity: O((V + E) log V), where V is the number of vertices and E is the number of edges.
 *   - Space Complexity: O(V + E), for MST storage.
 */

void TSPSolver::prim(Graph<int>* g, Vertex<int>* source, vector<Vertex<int>*> &result, double &cost) {
    Graph<int> * mst = new Graph<int>();
    MutablePriorityQueue<Vertex<int>> pq;
    for (auto v: g->getVertexSet()) {
        Vertex<int> * vMst = mst->findVertex(v->getInfo());
        if (!vMst) {
            mst->addVertex(v->getInfo(), 0.0, 0.0);
        }
        else {
            mst->addVertex(v->getInfo(), v->getLatitude(), v->getLatitude());
        }
        v->setVisited(false);
        v->setDist(INF);
        v->setPath(nullptr);
    }

    source->setDist(0);
    pq.insert(source);
    while (!pq.empty()) {
        Vertex<int>* u = pq.extractMin();
        if (u->isVisited()) {
            continue;
        }
        u->setVisited(true);
        if (u->getInfo() != source->getInfo()) {
            mst->addBidirectionalEdge(u->getPath()->getOrig()->getInfo(), u->getInfo(), u->getPath()->getWeight());
        }
        for (Edge<int>* e: u->getAdj()) {
            Vertex<int>* v = e->getDest();
            double w = e->getWeight();
            if (!v->isVisited() && w < v->getDist()) {
                double previous = v->getDist();
                v->setDist(w);
                v->setPath(e);
                if (previous == INF) {
                    pq.insert(v);
                } else {
                    pq.decreaseKey(v);
                }
            }
        }
    }
    for (auto v: g->getVertexSet()) {
        v->setVisited(false);
    }
    Vertex<int>* v = mst->findVertex(0);
    Vertex<int>* prev = nullptr;
    preorderMST(g, v, result, cost, prev);
    delete mst;
}

/**
 * @brief Calculates the TSP using a triangle inequality heuristic and prints the result.
 * @param g Pointer to the graph.
 * @details This function calculates the TSP path and cost using a triangle inequality heuristic (MST-based) and prints the results including elapsed time.
 *   - Time Complexity: O(V + E log V), where V is the number of vertices and E is the number of edges.
 *   - Space Complexity: O(V), for path storage.
 */

void TSPSolver::calculateTriangleTSP(Graph<int>* g) {
    vector<Vertex<int>*> tsp_path;
    double cost = 0;
    double costFinal = 0;
    auto start = std::chrono::high_resolution_clock::now();
    Vertex<int>* source = g->findVertex(0);


    prim(g, source, tsp_path, cost);
    for(auto e : tsp_path[tsp_path.size() -1]->getAdj()){
        if (e->getDest()->getInfo() == source->getInfo()){
            cost += e->getWeight();
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    Vertex<int>* finalVtx = g->findVertex(0);
    tsp_path.push_back(finalVtx);
    for (size_t i = 0; i < tsp_path.size()-1; i++) {
        auto nextNum = tsp_path[i+1]->getInfo();
        bool edgeExists = false;
        for(auto e : tsp_path[i]->getAdj()){
            if(e->getDest()->getInfo() == nextNum){
                costFinal += e->getWeight();
                edgeExists = true;
                break;
            }
        }
        if(!edgeExists && tsp_path[i]->getLatitude()){
            costFinal += haversineDistance(tsp_path[i], tsp_path[i+1]);
        }
    }

    cout << "TSP Path : ";
    for (size_t i = 0; i < tsp_path.size(); i++) {
        cout << tsp_path[i]->getInfo() << (i == tsp_path.size() - 1 ? "\n" : " -> ");
    }
    cout << '\n';

    cout << "Cost: " << costFinal << '\n';
    cout << "Elapsed Time: " << duration.count() << " s\n\n";
}

/**
 * @brief Performs a 2-opt swap on the TSP path between indices i and k.
 * @param tsp_path Reference to the vector storing the TSP path.
 * @param i The starting index for the swap.
 * @param k The ending index for the swap.
 * @details This function performs the 2-opt swap to potentially shorten the TSP path.
 *   - Time Complexity: O((k - i) / 2)
 *   - Space Complexity: O(1)
 */

void TSPSolver::twoOptSwap(vector<Vertex<int>*>& tsp_path, int i, int k) {
    while (i < k) {
        std::swap(tsp_path[i % tsp_path.size()], tsp_path[k % tsp_path.size()]);
        i++;
        k--;
    }
}

/**
 * @brief Improves the TSP path using the 2-opt algorithm.
 * @param tsp_path Reference to the vector storing the TSP path.
 * @param two_opt_iterations The number of iterations for the 2-opt algorithm.
 * @details This function iteratively applies the 2-opt swap to reduce the TSP path length.
 *   - Time Complexity: O(n^2 * iterations), where n is the number of vertices and iterations is the number of 2-opt iterations.
 *   - Space Complexity: O(1)
 */

void TSPSolver::twoOptAlgorithm(vector<Vertex<int>*>& tsp_path, unsigned int two_opt_iterations) {
    int n = tsp_path.size();
    unsigned int iterations = 0;
    bool improvement = true;

    while (improvement && iterations < two_opt_iterations) {
        iterations++;
        improvement = false;
        for (int i = 0; i < n - 2; ++i) {
            for (int k = i + 2; k < n; ++k) {
                auto a = tsp_path[i]->getEdge(tsp_path[i + 1]->getInfo());
                auto b = tsp_path[k]->getEdge(tsp_path[(k + 1) % n]->getInfo());
                auto c = tsp_path[i]->getEdge(tsp_path[k]->getInfo());
                auto d = tsp_path[i + 1]->getEdge(tsp_path[(k + 1) % n]->getInfo());
                if (a == nullptr || b == nullptr || c == nullptr || d == nullptr) {
                    continue;
                }
                double currentDistance = a->getWeight() + b->getWeight();
                double newDistance = c->getWeight() + d->getWeight();
                if (newDistance < currentDistance) {
                    twoOptSwap(tsp_path, i + 1, k);
                    improvement = true;
                }
            }
        }
    }
}

/**
 * @brief Solves the TSP using the nearest neighbor heuristic with optional 2-opt optimization.
 * @param g Pointer to the graph.
 * @param tsp_path Reference to the vector storing the TSP path.
 * @param two_opt_iterations The number of iterations for the 2-opt algorithm.
 * @return The total cost of the TSP path.
 * @details This function applies the nearest neighbor heuristic to generate an initial TSP path and then improves it using the 2-opt algorithm.
 *   - Time Complexity: O(V^2 + n^2 * iterations), where V is the number of vertices and iterations is the number of 2-opt iterations.
 *   - Space Complexity: O(V), for path storage.
 */

double TSPSolver::tspNearestNeighbor(Graph<int>* g, vector<Vertex<int>*>& tsp_path, unsigned int two_opt_iterations) {
    tsp_path.clear();
    std::size_t num_vertices = g->getNumVertex();
    std::vector<bool> visited(num_vertices, false);
    std::size_t start_idx = 0;
    Vertex<int>* current = g->getVertexSet()[start_idx];
    tsp_path.push_back(current);
    visited[start_idx] = true;

    while (tsp_path.size() < num_vertices) {
        double min_edge_weight = std::numeric_limits<double>::infinity();
        Vertex<int>* next_vertex = nullptr;

        for (Edge<int>* edge : current->getAdj()) {
            Vertex<int>* neighbor_vertex = edge->getDest();
            if (!visited[neighbor_vertex->getInfo()] && edge->getWeight() < min_edge_weight) {
                min_edge_weight = edge->getWeight();
                next_vertex = neighbor_vertex;
            }
        }

        if (next_vertex == nullptr) {
            break;
        }

        tsp_path.push_back(next_vertex);
        visited[next_vertex->getInfo()] = true;
        current = next_vertex;
    }

    if (tsp_path.size() == num_vertices) {
        Edge<int>* final_edge = tsp_path.back()->getEdge(g->getVertexSet()[start_idx]->getInfo());
        tsp_path.push_back(g->getVertexSet()[start_idx]);
    }
    twoOptAlgorithm(tsp_path, two_opt_iterations);

    double cost = 0;
    for (int i = 0; i < tsp_path.size() - 1; i++) {
        Edge<int>* edge = tsp_path[i]->getEdge(tsp_path[i + 1]->getInfo());
        cost += edge->getWeight();
    }

    return cost;
}

/**
 * @brief Calculates the TSP using the nearest neighbor heuristic and prints the result.
 * @param g Pointer to the graph.
 * @details This function calculates the TSP path and cost using the nearest neighbor heuristic with optional 2-opt optimization, and prints the results including elapsed time.
 *   - Time Complexity: O(V^2 + n^2 * iterations), where V is the number of vertices and iterations is the number of 2-opt iterations.
 *   - Space Complexity: O(V), for path storage.
 */

void TSPSolver::calculateNearestNeighborTSP(Graph<int>* g) {
    unsigned int iterations;
    cout << "Iterations for the 2opt algorithm : ";
    cin >> iterations;

    vector<Vertex<int>*> tsp_path;

    auto start = std::chrono::high_resolution_clock::now();

    double cost = tspNearestNeighbor(g,tsp_path, iterations);

    auto end = std::chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    cout << "TSP Path (Nearest Neighbor): ";
    for (int i = 0; i < tsp_path.size(); i++) {
        cout << tsp_path[i]->getInfo() << (i == tsp_path.size() - 1 ? "\n" : " -> ");
    }

    cout << "Cost: " << cost << '\n';
    cout << "Elapsed Time: " << duration.count() << " s\n\n";
}



/**
 * @brief Helper function for Prim's algorithm to find minimum spanning tree (MST).
 * @param g Pointer to the graph.
 * @param source Pointer to the source vertex.
 * @param result Reference to the vector to store the MST vertices.
 * @param cost Reference to the total cost of the MST.
 */
void TSPSolver::primHelper(Graph<int> *g, Vertex<int> *source, vector<Vertex<int> *> &result, double &cost) {
    Graph<int> *mst = new Graph<int>();
    MutablePriorityQueue<Vertex<int>> pq;

    // Initialize the MST graph and priority queue
    for (auto v : g->getVertexSet()) {
        Vertex<int> *vMst = mst->findVertex(v->getInfo());
        if (!vMst) {
            mst->addVertex(v->getInfo(), 0.0, 0.0);
        } else {
            mst->addVertex(v->getInfo(), v->getLatitude(), v->getLatitude());
        }
        v->setVisited(false);
        v->setDist(INF);
        v->setPath(nullptr);
    }

    // Initialize source vertex
    source->setDist(0);
    pq.insert(source);

    // Prim's algorithm
    while (!pq.empty()) {
        Vertex<int> *u = pq.extractMin();
        if (u->isVisited()) {
            continue;
        }
        u->setVisited(true);

        // Add edge to MST if not the source vertex
        if (u->getInfo() != source->getInfo()) {
            mst->addBidirectionalEdge(u->getPath()->getOrig()->getInfo(), u->getInfo(), u->getPath()->getWeight());
        }

        // Update distances and paths to adjacent vertices
        for (Edge<int> *e : u->getAdj()) {
            Vertex<int> *v = e->getDest();
            double w = e->getWeight();
            if (!v->isVisited() && w < v->getDist()) {
                double previous = v->getDist();
                v->setDist(w);
                v->setPath(e);
                if (previous == INF) {
                    pq.insert(v);
                } else {
                    pq.decreaseKey(v);
                }
            }
        }
    }

    // Reset visited status of vertices in the original graph
    for (auto v : g->getVertexSet()) {
        v->setVisited(false);
    }

    // Traverse the MST in preorder to construct the result
    Vertex<int> *v = g->getVertexSet()[0];
    Vertex<int> *prev = nullptr;
    preorderMST(g, v, result, cost, prev);

    delete mst;
}

/**
 * @brief Calculates the TSP path using Triangle TSP algorithm and returns the result.
 * @param g Pointer to the graph.
 * @return Vector of pointers to vertices representing the TSP path.
 */
vector<Vertex<int> *> TSPSolver::calculateTriangleTSPReturning(Graph<int> *g) {
    vector<Vertex<int> *> tsp_path;
    double cost = 0;

    // Start timer
    auto start = std::chrono::high_resolution_clock::now();
    Vertex<int> *source = g->getVertexSet()[0];

    // Run Prim's algorithm to find MST and construct TSP path
    primHelper(g, source, tsp_path, cost);

    // Adjust cost to account for returning to the source vertex
    for (auto e : tsp_path[tsp_path.size() - 1]->getAdj()) {
        if (e->getDest()->getInfo() == source->getInfo()) {
            cost += e->getWeight();
        }
    }

    // End timer
    auto end = std::chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    // Finalize the TSP path and calculate the total cost
    Vertex<int> *finalVtx = g->getVertexSet()[0];
    tsp_path.push_back(finalVtx);
    double costFinal = 0;
    for (size_t i = 0; i < tsp_path.size() - 1; i++) {
        auto nextNum = tsp_path[i + 1]->getInfo();
        bool edgeExists = false;
        for (auto e : tsp_path[i]->getAdj()) {
            if (e->getDest()->getInfo() == nextNum) {
                costFinal += e->getWeight();
                edgeExists = true;
                break;
            }
        }
        if (!edgeExists && tsp_path[i]->getLatitude()) {
            costFinal += haversineDistance(tsp_path[i], tsp_path[i + 1]);
        }
    }

    return tsp_path;
}


/**
 * @brief Initializes the cluster centers by selecting k random vertices from the graph.
 * @param clusters A reference to the vector of clusters to initialize.
 * @param graph A pointer to the graph containing the vertices.
 * @param k The number of clusters to initialize.
 */
void TSPSolver::initializeCenters(vector<Cluster>& clusters, const Graph<int>* graph, int k) {

    srand(time(NULL));
    vector<Vertex<int>*> vertices = graph->getVertexSet();
    int n = vertices.size();

    int index = rand() % n;
    const Vertex<int>* vertex = vertices[index];
    clusters[0].centerX = vertex->getLatitude();
    clusters[0].centerY = vertex->getLongitude();
    clusters[0].cities.push_back(vertex->getInfo());

    for (int i = 1; i < k; ++i) {
        vector<double> distances(n, numeric_limits<double>::max());

        for (int j = 0; j < n; ++j) {
            for (int l = 0; l < i; ++l) {
                double distance = haversineDistance(vertices[j], vertices[clusters[l].cities[0]]);
                distances[j] = min(distances[j], distance);
            }
        }

        double total = accumulate(distances.begin(), distances.end(), 0.0);
        double r = (double) rand() / RAND_MAX * total;
        double cumulative = 0.0;
        for (int j = 0; j < n; ++j) {
            cumulative += distances[j];
            if (cumulative >= r) {
                const Vertex<int>* selectedVertex = vertices[j];
                clusters[i].centerX = selectedVertex->getLatitude();
                clusters[i].centerY = selectedVertex->getLongitude();
                clusters[i].cities.push_back(selectedVertex->getInfo());
                break;
            }
        }
    }
}

/**
 * @brief Assigns each vertex in the graph to the nearest cluster center.
 * @param clusters A reference to the vector of clusters.
 * @param graph A pointer to the graph containing the vertices.
 */
void TSPSolver::assignToClusters(vector<Cluster>& clusters, const Graph<int>* graph) {
    for (Cluster& cluster : clusters)
        cluster.cities.clear();

    for (auto vertex : graph->getVertexSet()) {
        double minDistance = numeric_limits<double>::max();
        int nearestClusterIndex = -1;
        for (size_t i = 0; i < clusters.size(); ++i) {
            double distance = haversineDistance(vertex, graph->getVertexSet()[clusters[i].cities[0]]);
            if (distance < minDistance) {
                minDistance = distance;
                nearestClusterIndex = i;
            }
        }
        clusters[nearestClusterIndex].cities.push_back(vertex->getInfo());
    }
}

/**
 * @brief Updates the cluster centers based on the mean position of the vertices in each cluster.
 * @param clusters A reference to the vector of clusters.
 * @param graph A pointer to the graph containing the vertices.
 */
void TSPSolver::updateCenters(vector<Cluster>& clusters, const Graph<int>* graph) {
    for (Cluster& cluster : clusters) {
        double sumX = 0.0;
        double sumY = 0.0;
        for (int cityIndex : cluster.cities) {
            const Vertex<int>* vertex = graph->findVertex(cityIndex);
            sumX += vertex->getLatitude();
            sumY += vertex->getLongitude();
        }
        cluster.centerX = sumX / cluster.cities.size();
        cluster.centerY = sumY / cluster.cities.size();
    }
}

/**
 * @brief Performs k-means clustering on the graph vertices.
 * @param graph A pointer to the graph containing the vertices.
 * @param k The number of clusters to form.
 * @param maxIterations The maximum number of iterations for the k-means algorithm.
 * @return A vector of clusters after the k-means algorithm has been applied.
 */
vector<Cluster> TSPSolver::kMeansClustering(const Graph<int>* graph) {

    printClustersHelper();

    int k, maxIterations;
    cout << "Insert number of clusters: ";
    cin >> k;
    cout << endl;
    cout << "Insert number of iterations: ";
    cin >> maxIterations;
    this_thread::sleep_for(chrono::milliseconds(1000));
    cout << endl;
    cout << "This will only take a few moments..." << endl << endl;

    //k = 100;
    //maxIterations = 1000;
    std::vector<Cluster> clusters(k);

    initializeCenters(clusters, graph, k);

    for (int iter = 0; iter < maxIterations; ++iter) {

        assignToClusters(clusters, graph);

        updateCenters(clusters, graph);
    }
    return clusters;
}

/**
 * @brief Finds the best tour for a given cluster using the TSP algorithm.
 * @param graph A pointer to the original graph containing the vertices.
 * @param cluster A reference to the cluster for which to find the best tour.
 * @return A vector of pointers to vertices representing the best tour for the cluster.
 */


vector<Vertex<int>*> TSPSolver::findBestTourForCluster(const Graph<int> *graph, const Cluster &cluster, bool two) {
    Graph<int> subgraph;
    unordered_map<int, Vertex<int>*> vertexMap;

    for (int id : cluster.cities) {
        auto vertex = graph->findVertex(id);
        subgraph.addVertex(id, vertex->getLongitude(), vertex->getLatitude());
        vertexMap[id] = vertex;
    }
    for (int id : cluster.cities) {
        auto vertex = graph->findVertex(id);
        for (const auto& edge : vertex->getAdj()) {
            int destCityId = edge->getDest()->getInfo();
            if (vertexMap.find(destCityId) != vertexMap.end())
                subgraph.addEdge(vertex->getInfo(), destCityId, edge->getWeight());
        }
    }

    return two ? optimizeTourTwoOpt(calculateTriangleTSPReturning(&subgraph)) :  calculateTriangleTSPReturning(&subgraph);

}

/**
 * @brief Finds the closest city to a given city from a tour.
 * @param city A pointer to the city for which to find the closest city.
 * @param tour A vector of pointers to vertices representing the tour.
 * @return A pointer to the closest city in the tour.
 */
Vertex<int>* TSPSolver::findClosestCity(const Graph<int>* graph, Vertex<int> *city, const std::vector<Vertex<int>*> &tour) {
    Vertex<int>* closestCity = nullptr;
    double minDistance = std::numeric_limits<double>::max();

    for (Vertex<int>* v : tour) {
        double distance = haversineDistance(city, v);
        if (distance < minDistance) {
            minDistance = distance;
            closestCity = v;
        }
    }

    return closestCity;
}

/**
 * @brief Unites all cluster tours into a single complete tour.
 * @param graph A pointer to the original graph containing the vertices.
 * @param clusters A reference to the vector of clusters.
 * @return A vector of pointers to vertices representing the complete tour.
 */
void TSPSolver::uniteAllClusterTours(const Graph<int>* graph, const vector<Cluster>& clusters, bool two) {

    auto start = std::chrono::high_resolution_clock::now();
    vector<Vertex<int>*> completeTour;

    vector<std::vector<Vertex<int>*>> clusterTours;
    for (const auto& cluster : clusters) {
        std::vector<Vertex<int>*> tour = findBestTourForCluster(graph, cluster, two);
        clusterTours.push_back(tour);
    }

    completeTour = clusterTours[0];

    for (size_t i = 1; i < clusterTours.size(); ++i) {
        Vertex<int>* closestCity = findClosestCity(graph,completeTour.back(), clusterTours[i]);
        auto it = find(clusterTours[i].begin(), clusterTours[i].end(), closestCity);
        completeTour.insert(completeTour.end(), std::next(it), clusterTours[i].end());
        completeTour.insert(completeTour.end(), clusterTours[i].begin(), it);
    }

    completeTour.push_back(completeTour.front());

    auto end = std::chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    double totalCost = 0.0;
    cout << "TSP Path : ";
    for (size_t i = 0; i < completeTour.size() - 1; i++) {
        cout << completeTour[i]->getInfo() << (i == completeTour.size() - 1 ? "\n" : " -> ");
        for (auto edge : completeTour[i]->getAdj()) {
            if (edge->getDest()->getInfo() == completeTour[i + 1]->getInfo())
                totalCost += edge->getWeight();
        }
    }
    cout << completeTour[completeTour.size() - 1]->getInfo() << endl << endl;
    cout << "Total cost: " << totalCost << endl;
    cout << '\n';

    cout << "Elapsed Time: " << duration.count() << " s\n\n";

}


vector<Vertex<int>*> TSPSolver::twoOptSwapHelper(const vector<Vertex<int>*>& tour, int i, int k) {
    vector<Vertex<int>*> newTour(tour.size());
    // 1. take route[0] to route[i-1] and add them in order to new_route
    for (int c = 0; c <= i - 1; ++c) {
        newTour[c] = tour[c];
    }
    // 2. take route[i] to route[k] and add them in reverse order to new_route
    for (int c = i; c <= k; ++c) {
        newTour[c] = tour[k - (c - i)];
    }
    // 3. take route[k+1] to end and add them in order to new_route
    for (int c = k + 1; c < tour.size(); ++c) {
        newTour[c] = tour[c];
    }
    return newTour;
}

vector<Vertex<int>*> TSPSolver::optimizeTourTwoOpt(const vector<Vertex<int>*>& tour) {
    vector<Vertex<int>*> bestTour = tour;
    double bestDistance = calculateTourDistance(tour);
    bool improvement = true;

    while (improvement) {
        improvement = false;
        for (int i = 1; i < tour.size() - 1; ++i) {
            for (int k = i + 1; k < tour.size() - 1; ++k) {
                vector<Vertex<int>*> newTour = twoOptSwapHelper(bestTour, i, k);
                double newDistance = calculateTourDistance(newTour);
                if (newDistance < bestDistance) {
                    bestTour = newTour;
                    bestDistance = newDistance;
                    improvement = true;
                }
            }
        }
    }
    return bestTour;
}


double TSPSolver::calculateTourDistance(const vector<Vertex<int>*>& tour) {
    double totalDistance = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        totalDistance += haversineDistance(tour[i], tour[i + 1]);
    }
    return totalDistance;
}


