#include <iostream>
#include <iomanip>

#include "headers/Script.h"
#include "headers/TSPSolver.h"
#include "headers/Menu.h"

using namespace std;

int main() {
    Script script;
    TSPSolver tspSolver;
    Menu menu = Menu(tspSolver, script);
    menu.chooseGraph();
    return 0;
}
