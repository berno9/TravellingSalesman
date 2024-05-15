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
                if (edg->getDest()->getInfo() == e->getDest()->getInfo())
                    tmp = edg->getWeight();
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
        mst->addVertex(v->getInfo(), v->getLatitude(), v->getLatitude());
        v->setVisited(false);
        v->setDist(std::numeric_limits<double>::max());
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
            mst->addBidirectionalEdge(u->getPath()->getOrig()->getInfo(), u->getInfo(), u->getPath()->getWeight());//////////////
        }
        for (Edge<int>* e: u->getAdj()) {
            Vertex<int>* v = e->getDest();
            double w = e->getWeight();
            if (!v->isVisited() && w < v->getDist()) {
                double previous = v->getDist();
                v->setDist(w);
                v->setPath(e);
                if (previous == std::numeric_limits<double>::max()) {
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
        for(auto e : tsp_path[i]->getAdj()){
            if(e->getDest()->getInfo() == nextNum){
                cost2 += e->getWeight();
            }
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





