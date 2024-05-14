#include "../headers/TSPSolver.h"


Graph<int> *TSPSolver::TheG() {
    return script.getRealWorldGraph1();
}

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

    // Haversine formula
    double a = sin(dlat / 2) * sin(dlat / 2) +
               cos(lat1) * cos(lat2) *
               sin(dlon / 2) * sin(dlon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    double distance = R * c;

    return distance;
}


void TSPSolver::dijkstra(Graph<int>* g, Vertex<int>* source) {
    auto cmp = [](Vertex<int>* a, Vertex<int>* b) {
        return a->getDist() < b->getDist();
    };
    priority_queue<Vertex<int> *, std::vector<Vertex<int> *>, decltype(cmp)> pq(cmp);

    for (auto v: g->getVertexSet()) {
        v->setVisited(false);
        v->setDist(INF);
    }

    source->setDist(0);
    pq.push(source);
    while (!pq.empty()) {
        Vertex<int>* u = pq.top(); pq.pop();
        u->setVisited(true);

        for (Edge<int>* e: u->getAdj()) {
            Vertex<int>* v = e->getDest();
            double w = e->getWeight();
            if (!v->isVisited() && u->getDist() != INF && u->getDist() + w < v->getDist()) {
                v->setDist(u->getDist() + w);
                pq.push(v);
            }
        }
    }
}

void TSPSolver::tspBacktrack(Graph<int>* g, Vertex<int> *current, double current_cost, int num_visited, double &min_cost,
                                       vector<Vertex<int> *> &tsp_path) {
    if (num_visited == g->getNumVertex()) {
        double cost = current_cost;
        bool hasEdge = false;
        for (Edge<int> *e: current->getAdj()) {
            if (e->getDest()->getInfo() == 0) {
                hasEdge = true;
                cost += e->getWeight();
                break;
            }
        }

        if (!hasEdge) {
            dijkstra(g,current);
            cost += g->findVertex(current->getInfo())->getDist();
        }
        if (cost < min_cost) {
            min_cost = cost;

            Vertex<int> *init = g->findVertex(0);

            tsp_path.clear();
            tsp_path.push_back(init);
            tsp_path.push_back(current);
            for (Edge<int> *e = current->getPath();
                 e->getOrig()->getInfo() != g->findVertex(0)->getInfo(); e = e->getOrig()->getPath()) {
                tsp_path.push_back(e->getOrig());
            }
            tsp_path.push_back(init);

            reverse(tsp_path.begin(), tsp_path.end());
        }
        return;
    }

    for (Edge<int>* e: current->getAdj()) {
        Vertex<int>* w = e->getDest();
        if (!w->isVisited()) {
            w->setVisited(true);
            w->setPath(e);
            tspBacktrack(g,w, current_cost + e->getWeight(), num_visited + 1, min_cost, tsp_path);
            w->setVisited(false);
            w->setPath(nullptr);
        }
    }

}

double TSPSolver::tspBruteForce(Graph<int>* g, vector<Vertex<int> *> &tsp_path) {
    for (auto v: g->getVertexSet()) {
        v->setVisited(false);
        v->setPath(nullptr);
    }

    double min_cost = INF;
    Vertex<int>* source = g->findVertex(0);
    source->setVisited(true);
    tspBacktrack(g,source, 0, 1, min_cost, tsp_path);

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





void TSPSolver::initializeCenters(std::vector<Cluster>& clusters, const Graph<int>* graph, int k) {
    srand(time(NULL));
    for (int i = 0; i < k; ++i) {
        int index = rand() % graph->getVertexSet().size();
        const Vertex<int>* vertex = graph->getVertexSet()[index];
        clusters[i].centerX = vertex->getLatitude();
        clusters[i].centerY = vertex->getLongitude();
    }
}

void TSPSolver::assignToClusters(std::vector<Cluster>& clusters, const Graph<int>* graph) {
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

void TSPSolver::updateCenters(std::vector<Cluster>& clusters, const Graph<int>* graph) {
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

std::vector<Cluster> TSPSolver::kMeansClustering(const Graph<int>* graph, int k, int maxIterations) {
    std::vector<Cluster> clusters(k);

    initializeCenters(clusters, graph, k);

    // Iterate until convergence or maximum iterations reached
    for (int iter = 0; iter < maxIterations; ++iter) {

        assignToClusters(clusters, graph);

        updateCenters(clusters, graph);

        for (Cluster& cluster : clusters) {
            cluster.cities.clear();
        }
    }

    return clusters;
}




