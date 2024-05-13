
#ifndef TRAVELLINGSALESMAN_TSPSOLVER_H
#define TRAVELLINGSALESMAN_TSPSOLVER_H

#include "Graph.h"
#include <vector>
#include <limits>
#include <functional>
#include <chrono>
#include <unordered_set>
#include <cmath>

using namespace std;

class TSPSolver{
private:


public:

    void dijkstra(Graph<int>* g, Vertex<int>* source);
    void tspBacktrack(Graph<int>* g, Vertex<int>* current, double current_cost, int numVisited,
                      double& minCost, vector<Vertex<int> *> &tspPath);
    double tspBruteForce(Graph<int>* g, vector<Vertex<int> *> &tspPath);
    void calculateTSP(Graph<int>*g);

    void triangularApproximation(Graph<int>* g, Vertex<int>* current, double current_cost, int num_visited, double& min_cost,
                                            vector<Vertex<int>*>& tsp_path);
    double tspTriangleBruteForce(Graph<int>* g, vector<Vertex<int> *> &tsp_path);
    void calculateTriangleTSP(Graph<int>* g) ;
    double haversineDistance(Vertex<int> *v1, Vertex<int> *v2);

    void prim(Graph<int>* g,Vertex<int>* source, vector<Vertex<int>*> &result, double &cost);
    void preorderMST(Graph<int>* g, Vertex<int>* current, std::vector<Vertex<int>*> &result, double &cost, Vertex<int>* &prev);

};


#endif //TRAVELLINGSALESMAN_TSPSOLVER_H
