
#ifndef TRAVELLINGSALESMAN_TSPSOLVER_H
#define TRAVELLINGSALESMAN_TSPSOLVER_H

#include "Graph.h"
#include <vector>
#include <limits>
#include <functional>
#include <chrono>
#include <unordered_set>

using namespace std;

class TSPSolver{
private:


public:

    void dijkstra(Graph<int>* g, Vertex<int>* source);
    void tspBacktrack(Graph<int>* g, Vertex<int>* current, double current_cost, int numVisited,
                      double& minCost, vector<Vertex<int> *> &tspPath);
    double tspBruteForce(Graph<int>* g, vector<Vertex<int> *> &tspPath);
    void calculateTSP(Graph<int>*g);
};


#endif //TRAVELLINGSALESMAN_TSPSOLVER_H
