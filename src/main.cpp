#include <iostream>
#include <iomanip>

#include "headers/Script.h"
#include "headers/TSPSolver.h"
#include "headers/Menu.h"

using namespace std;

int main() {
    Script script;
    TSPSolver tspSolver;
    //Menu menu = Menu(tspSolver, script);
    //menu.chooseGraph();

    // ------Testing------
    script.read_rwg("../cmake-build-debug/datasets/Real-World Graphs/graph1/nodes.csv", "../cmake-build-debug/datasets/Real-World Graphs/graph1/edges.csv");
    auto g = script.getRealWorldGraph();
    tspSolver.uniteAllClusterTours(g, tspSolver.kMeansClustering(g));


    return 0;
}
