
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


    double haversineDistance(Vertex<int> *v1, Vertex<int> *v2);

    void prim(Graph<int>* g,Vertex<int>* source, vector<Vertex<int>*> &result, double &cost);
    void preorderMST(Graph<int>* g, Vertex<int>* current, std::vector<Vertex<int>*> &result, double &cost, Vertex<int>* &prev);


    void tspBacktrack(Graph<int>* g, Vertex<int>* current, double current_cost, int num_visited, double &min_cost, vector<Vertex<int> *> &tsp_path, vector<Vertex<int> *> &current_path);
    double tspBruteForce(Graph<int>* g, vector<Vertex<int> *> &tspPath);
    void calculateTSP(Graph<int>*g);

    void triangularApproximation(Graph<int>* g, Vertex<int>* current, double current_cost, int num_visited, double& min_cost,
                                 vector<Vertex<int>*>& tsp_path);
    double tspTriangleBruteForce(Graph<int>* g, vector<Vertex<int> *> &tsp_path);
    void calculateTriangleTSP(Graph<int>* g);
    vector<Vertex<int>*> calculateTriangleTSPReturning(Graph<int>* g);


    void initializeCenters(vector<Cluster>& clusters, const Graph<int>* graph, int k); // Initialize the cluster centers randomly
    void assignToClusters(vector<Cluster>& clusters, const Graph<int>* graph); // Assign vertices (cities) to clusters based on the nearest cluster center
    void updateCenters(vector<Cluster>& clusters, const Graph<int>* graph); // Update the cluster centers based on the mean of the vertices in each cluster
    vector<Cluster> kMeansClustering(const Graph<int>* graph, int k, int maxIterations); // Perform KMeans clustering

    vector<Vertex<int>*> findBestTourForCluster(const Graph<int>* graph, const Cluster& cluster); // Find best tour on a single cluster
    Vertex<int>* findClosestCity(Vertex<int>* city, const vector<Vertex<int>*>& tour);  // Find closest city to a given tour
    vector<Vertex<int>*> uniteAllClusterTours(const Graph<int>* graph, const vector<Cluster>& clusters);  // Find complete tour from all clusters


    void twoOptSwap(vector<Vertex<int>*>& tsp_path, int i, int k);
    void twoOptAlgorithm(vector<Vertex<int>*>& tsp_path, unsigned int two_opt_iterations);
    double tspNearestNeighbor(Graph<int>* g, vector<Vertex<int>*>& tsp_path, unsigned int two_opt_iterations);
    void calculateNearestNeighborTSP(Graph<int>* g);


};


#endif //TRAVELLINGSALESMAN_TSPSOLVER_H
