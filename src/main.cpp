#include <iostream>

#include "headers/Script.h"
#include "headers/TSPSolver.h"
#include "headers/Menu.h"

using namespace std;

int main() {
    Script script;
    script.read_tourism();
    //script.read_shipping();
    //script.read_tourism();
    TSPSolver tspSolver;
    tspSolver.calculateTSP(script.getTmGraph());

    Menu menu = Menu(tspSolver, script);
    //menu.mainMenu();
    return 0;
}
