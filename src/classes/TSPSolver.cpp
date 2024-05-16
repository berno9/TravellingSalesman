#include "../headers/TSPSolver.h"


Graph<int> *TSPSolver::TheG() {
    return script.getRealWorldGraph1();
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
    cout << "Iterations for the 2opt algorithm (0 to not run 2opt algorithm, more iterations -> more precise): ";
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




/*


void TSPSolver::initializeCenters(vector<Cluster>& clusters, const Graph<int>* graph, int k) {
    srand(time(NULL));
    unordered_map<int, bool> usedIndices;

    for (int i = 0; i < k; ++i) {
        int index;
        do {
            index = rand() % graph->getVertexSet().size();
        } while (usedIndices.find(index) != usedIndices.end());

        usedIndices[index] = true;
        const Vertex<int>* vertex = graph->getVertexSet()[index];
        clusters[i].centerX = vertex->getLatitude();
        clusters[i].centerY = vertex->getLongitude();
        clusters[i].cities.push_back(vertex->getInfo());
    }
}

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

vector<Cluster> TSPSolver::kMeansClustering(const Graph<int>* graph, int k, int maxIterations) {
    std::vector<Cluster> clusters(k);

    initializeCenters(clusters, graph, k);

    // Iterate until convergence or maximum iterations reached
    for (int iter = 0; iter < maxIterations; ++iter) {

        assignToClusters(clusters, graph);

        updateCenters(clusters, graph);
    }
    return clusters;
}

vector<Vertex<int> *> TSPSolver::findBestTourForCluster(const Graph<int> *graph, const Cluster &cluster) {

    Graph<int> subgraph;
    std::unordered_map<int, Vertex<int>*> vertexMap;

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

    return calculateTriangleTSPReturning(&subgraph);
}

Vertex<int> *TSPSolver::findClosestCity(Vertex<int> *city, const std::vector<Vertex<int> *> &tour) {
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

vector<Vertex<int>*> TSPSolver::uniteAllClusterTours(const Graph<int>* graph, const vector<Cluster>& clusters) {
    vector<Vertex<int>*> completeTour;

    // Step 1: Calculate TSP Tour for Each Cluster
    vector<std::vector<Vertex<int>*>> clusterTours;
    for (const auto& cluster : clusters) {
        std::vector<Vertex<int>*> tour = findBestTourForCluster(graph, cluster);
        clusterTours.push_back(tour);
    }

    // Step 2: Tour stitching
    completeTour = clusterTours[0];

    // Stitch the TSP tours of the remaining clusters
    for (size_t i = 1; i < clusterTours.size(); ++i) {
        Vertex<int>* closestCity = findClosestCity(completeTour.back(), clusterTours[i]);
        auto it = find(clusterTours[i].begin(), clusterTours[i].end(), closestCity);
        completeTour.insert(completeTour.end(), std::next(it), clusterTours[i].end());
        completeTour.insert(completeTour.end(), clusterTours[i].begin(), it);
    }

    // Step 3: Close the Tour
    completeTour.push_back(completeTour.front());

    return completeTour;
}

*/
