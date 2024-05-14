#include <iostream>
#include <iomanip>

#include "headers/Script.h"
#include "headers/TSPSolver.h"
#include "headers/Menu.h"

using namespace std;

int main() {
    Script script;
    //script.read_tourism();
    //script.read_shipping();
    //script.read_tourism();
    TSPSolver tspSolver;
    //tspSolver.calculateTSP(script.getTmGraph());

    //script.read_stadiums();
    //script.read_shipping();
    //script.read_tourism();
    script.read_rwg_g1();
    //TSPSolver tspSolver;
    //tspSolver.calculateTSP(script.getStGraph());


    //Menu menu = Menu(tspSolver, script);
    //menu.mainMenu();

    auto g1 = script.getRealWorldGraph1();
    cout << "ola" << endl;
    int k = 4; // number of expected clusters
    cout << "ola" << endl;
    int maxIterations = 70; // max number of iterations to achieve convergence
    cout << "ola" << endl;
    auto sol = tspSolver.kMeansClustering(g1, k, maxIterations);
    cout << "ola" << endl;

    for (auto s : sol) {
        cout << "One cluster: " << s.centerX << ", " << s.centerY << endl;
        for (auto e : s.cities)
            cout << e << endl;
        cout << "Next cluster!" << endl;
    }
    cout << "ola" << endl;

    /*for (auto v : g1->getVertexSet()) {
        cout << v->getInfo() << ", " << fixed << setprecision(6) << v->getLongitude() << ", "
                << fixed << setprecision(6) << v->getLatitude() << "; " << endl;
    }*/


    return 0;
}
