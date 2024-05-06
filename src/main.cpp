#include <iostream>
#include <iomanip>

#include "headers/Script.h"
#include "headers/TSPSolver.h"
#include "headers/Menu.h"

using namespace std;

int main() {
    Script script;
    //script.read_stadiums();
    //script.read_shipping();
    //script.read_tourism();
    script.read_rwg_g1();
    //TSPSolver tspSolver;
    //tspSolver.calculateTSP(script.getStGraph());

    //Menu menu = Menu(tspSolver, script);
    //menu.mainMenu();

    auto g1 = script.getRealWorldGraph1();

    for (auto v : g1->getVertexSet()) {
        cout << v->getInfo() << ", " << fixed << setprecision(6) << v->getLongitude() << ", "
                << fixed << setprecision(6) << v->getLatitude() << "; " << endl;
    }
    return 0;
}
