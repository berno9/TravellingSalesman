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



void TSPSolver::triangularApproximation(Graph<int>* g, Vertex<int>* current, double current_cost, int num_visited, double& min_cost,
                             vector<Vertex<int>*>& tsp_path) {
    if (num_visited == g->getNumVertex()) {
        double cost = current_cost;
        bool hasEdge = false;
        for (Edge<int>* e : current->getAdj()) {
            if (e->getDest()->getInfo() == 0) {
                hasEdge = true;
                cost += e->getWeight();
                break;
            }
        }

        if (!hasEdge) {
            // Find the nearest unvisited vertex that connects to the starting vertex
            Vertex<int>* nearestStart = nullptr;
            double minStartDist = INF;
            for (auto v : g->getVertexSet()) {
                if (!v->isVisited() && v->getDist() < minStartDist) {
                    nearestStart = v;
                    minStartDist = v->getDist();
                }
            }

            // Find the nearest unvisited vertex that connects to the current vertex
            Vertex<int>* nearestCurrent = nullptr;
            double minCurrentDist = INF;
            for (Edge<int>* e : current->getAdj()) {
                Vertex<int>* neighbor = e->getDest();
                if (!neighbor->isVisited() && neighbor->getDist() < minCurrentDist) {
                    nearestCurrent = neighbor;
                    minCurrentDist = neighbor->getDist();
                }
            }

            // If there's a suitable intermediate vertex, update the cost
            if (nearestStart && nearestCurrent &&
                current->getDist() + nearestCurrent->getDist() == nearestCurrent->getDist()) {
                dijkstra(g, current); // Update distances from the current vertex
                cost += nearestCurrent->getDist(); // Add the distance from the intermediate vertex to the starting vertex
            } else {
                dijkstra(g, current); // Update distances from the current vertex
                cost += g->findVertex(current->getInfo())->getDist(); // Add the distance to the starting vertex
            }
        }

        if (cost < min_cost) {
            min_cost = cost;

            Vertex<int>* init = g->findVertex(0);

            tsp_path.clear();
            tsp_path.push_back(init);
            tsp_path.push_back(current);
            for (Edge<int>* e = current->getPath(); e->getOrig()->getInfo() != g->findVertex(0)->getInfo();
                 e = e->getOrig()->getPath()) {
                tsp_path.push_back(e->getOrig());
            }
            tsp_path.push_back(init);

            reverse(tsp_path.begin(), tsp_path.end());
        }
        return;
    }

    for (Edge<int>* e : current->getAdj()) {
        Vertex<int>* w = e->getDest();
        if (!w->isVisited()) {
            w->setVisited(true);
            w->setPath(e);
            triangularApproximation(g, w, current_cost + e->getWeight(), num_visited + 1, min_cost, tsp_path);
            w->setVisited(false);
            w->setPath(nullptr);
        }
    }
}

double TSPSolver::tspTriangleBruteForce(Graph<int>* g, vector<Vertex<int> *> &tsp_path) {
    for (auto v: g->getVertexSet()) {
        v->setVisited(false);
        v->setPath(nullptr);
    }

    double min_cost = INF;
    Vertex<int>* source = g->findVertex(0);
    source->setVisited(true);
    triangularApproximation(g,source, 0, 1, min_cost, tsp_path);

    return min_cost;
}
void TSPSolver::calculateTriangleTSP(Graph<int>* g) {
    vector<Vertex<int> *> tsp_path;

    auto start = std::chrono::high_resolution_clock::now();
    double cost = tspTriangleBruteForce(g, tsp_path);
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





