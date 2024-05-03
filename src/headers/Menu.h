
#ifndef TRAVELLINGSALESMAN_MENU_H
#define TRAVELLINGSALESMAN_MENU_H


#include "TSPSolver.h"
#include "Script.h"

class Menu {
private:
    Script script;
    TSPSolver tspSolver;
public:
    Menu(TSPSolver tspSolver, Script script);

    void mainMenu();
    void switchGraph();
    void printToyGraph();

    void solutionTSPBack(Graph<int>* g);
    void next();
};


#endif //TRAVELLINGSALESMAN_MENU_H
