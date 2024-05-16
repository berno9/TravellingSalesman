#include <iostream>
#include <iomanip>

#include "headers/Script.h"
#include "headers/TSPSolver.h"
#include "headers/Menu.h"

using namespace std;

int main() {
    Script script;
    TSPSolver tspSolver;
    script.read_rwg_g1();

    /*
    script.read_stadiums();
    script.read_shipping();
    script.read_tourism();
    tspSolver.calculateTSP(script.getTmGraph());
    tspSolver.calculateTriangleTSP(script.getTmGraph());

     */
    //Menu menu = Menu(tspSolver, script);
    //menu.chooseGraph();
    //tspSolver.calculateTSP(script.getStGraph());
    tspSolver.calculateTriangleTSP(script.getRealWorldGraph1());
    //tspSolver.calculateNearestNeighborTSP(script.getRealWorldGraph1());


//---------------------Bernardo----------------------------
    /*
    auto g1 = script.getRealWorldGraph1();

    for (auto v : g1->getVertexSet()) {
        cout << v->getInfo() << ", " << fixed << setprecision(6) << v->getLongitude() << ", "
                << fixed << setprecision(6) << v->getLatitude() << "; " << endl;
    }*/
//--------------------------------------------------------
    return 0;
}
