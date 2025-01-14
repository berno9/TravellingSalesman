

#ifndef TRAVELLINGSALESMAN_MENU_H
#define TRAVELLINGSALESMAN_MENU_H


#include "TSPSolver.h"
#include "Script.h"

class Menu {
private:
    Graph<int> *g;
    Script script;
    TSPSolver tspSolver;
public:
    Menu(TSPSolver tspSolver, Script script);
    void mainMenu();
    void chooseGraph();
    void printToyGraph();
    void printRealWorldGraph();
    void printExtraFullyConnected();

    void ifTwoOptHeuristic();

};


#endif //TRAVELLINGSALESMAN_MENU_H

