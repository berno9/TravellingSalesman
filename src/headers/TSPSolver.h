
#ifndef TRAVELLINGSALESMAN_TSPSOLVER_H
#define TRAVELLINGSALESMAN_TSPSOLVER_H

#include "Graph.h"
#include <vector>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <unordered_set>
#include "Script.h"


using namespace std;


struct Cluster {
    double centerX;
    double centerY;
    std::vector<int> cities;
};

class TSPSolver{
private:
    Script script;

public:

    Graph<int>* TheG();

    void dijkstra(Graph<int>* g, Vertex<int>* source);
    void tspBacktrack(Graph<int>* g, Vertex<int>* current, double current_cost, int numVisited,
                      double& minCost, vector<Vertex<int> *> &tspPath);
    double tspBruteForce(Graph<int>* g, vector<Vertex<int> *> &tspPath);
    void calculateTSP(Graph<int>*g);

    double haversineDistance(Vertex<int> *v1, Vertex<int> *v2);


    // alternative TSP heuristic
    // first, implement clusters
    void initializeCenters(std::vector<Cluster>& clusters, const Graph<int>* graph, int k); // Function to initialize the cluster centers randomly
    void assignToClusters(std::vector<Cluster>& clusters, const Graph<int>* graph); // Function to assign vertices (cities) to clusters based on the nearest cluster center
    void updateCenters(std::vector<Cluster>& clusters, const Graph<int>* graph); // Function to update the cluster centers based on the mean of the vertices in each cluster
    std::vector<Cluster> kMeansClustering(const Graph<int>* graph, int k, int maxIterations); // Function to perform KMeans clustering

};


#endif //TRAVELLINGSALESMAN_TSPSOLVER_H
