
#include "../headers/TSPSolver.h"

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


