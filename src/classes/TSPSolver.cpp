#include "../headers/TSPSolver.h"


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


void TSPSolver::prim(Graph<int>* g,Vertex<int>* source, vector<Vertex<int>*> &result, double &cost) {
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


void TSPSolver::calculateTriangleTSP(Graph<int>* g) {
    vector<Vertex<int>*> tsp_path;
    double cost = 0;
    double cost2 = 0;
    auto start = std::chrono::high_resolution_clock::now();
    Vertex<int>* source = g->findVertex(0);
    bool real_world = false;


    prim(g,source, tsp_path, cost);
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
                cost2 += e->getWeight();
                edgeExists = true;
                break;
            }
        }
        if(!edgeExists && tsp_path[i]->getLatitude()){
            cost2 += haversineDistance(tsp_path[i], tsp_path[i+1]);
        }
    }

    cout << "TSP Path : ";
    for (size_t i = 0; i < tsp_path.size(); i++) {
        cout << tsp_path[i]->getInfo() << (i == tsp_path.size() - 1 ? "\n" : " -> ");
    }
    cout << '\n';
    cout << "Costi: " << cost << '\n';
    cout << "Cost2: " << cost2 << '\n';
    cout << "Elapsed Time: " << duration.count() << " s\n\n";
}

void TSPSolver::twoOptSwap(vector<Vertex<int>*>& tsp_path, int i, int k) {
    while (i < k) {
        std::swap(tsp_path[i % tsp_path.size()], tsp_path[k % tsp_path.size()]);
        i++;
        k--;
    }
}

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
