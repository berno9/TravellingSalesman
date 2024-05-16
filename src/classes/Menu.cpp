

#include "../headers/Menu.h"

Menu::Menu(TSPSolver tspSolver, Script script) : tspSolver(tspSolver), script(script) {}

void Menu::mainMenu() {
    std::cout << std::endl;
    std::cout << "########################################################################################" << std::endl;
    std::cout << "##                                                                                    ##" << std::endl;
    std::cout << "##                    * * *   Consult TSP solutions   * * *                           ##" << std::endl;
    std::cout << "##                                                                                    ##" << std::endl;
    std::cout << "##      1 - Backtracking                                                              ##" << std::endl;
    std::cout << "##                                                                                    ##" << std::endl;
    std::cout << "##      2 - Triangular Approximation Heuristic                                        ##" << std::endl;
    std::cout << "##                                                                                    ##" << std::endl;
    std::cout << "##      3 - Divide cities into clusters and merge all sub tours                       ##" << std::endl;
    std::cout << "##                                                                                    ##" << std::endl;
    std::cout << "##      4 - Nearest neighbours                                                        ##" << std::endl;
    std::cout << "##                                                                                    ##" << std::endl;
    std::cout << "##      5 - Switch current graph                                                      ##" << std::endl;
    std::cout << "##                                                                                    ##" << std::endl;
    std::cout << "##      0-> Exit                                                                      ##" << std::endl;
    std::cout << "##                                                                                    ##" << std::endl;
    std::cout << "########################################################################################" << std::endl << std::endl;


    int k;
    std::cout << "Option: ";
    std::cin >> k;

    if (cin.fail()) {
        cin.clear();
        cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        cout << "Invalid input. Please enter an integer." << std::endl;
        mainMenu();
        return;
    }

    switch (k) {
        case 0:
            exit(0);
        case 1:
            // backtracking
            tspSolver.calculateTSP(g);
            mainMenu();
            break;
        case 2:
            // 2-approximation
            tspSolver.calculateTriangleTSP(g);
            mainMenu();
            break;
        case 3:
            // clusters
            mainMenu();
            break;
        case 4:
            // nearestNeighbor
            tspSolver.calculateNearestNeighborTSP(g);
            mainMenu();
            break;
        case 5:
            // nearestNeighbor
            chooseGraph();
            break;
        default:
            std::cout << "Please insert a valid option." << std::endl;
            mainMenu(); break;

    }
}

void Menu::chooseGraph() {
    std::cout << std::endl;
    std::cout << "#######################################################################" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##                 * * *   Type of Graph   * * *                     ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      1 - Toy-Graphs                                               ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      2 - Extra Fully Connected Graphs                             ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      2 - Real World Graphs                                        ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      0-> Exit                                                     ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "#######################################################################" << std::endl << std::endl;


    int k;
    std::cout << "Option: ";
    std::cin >> k;

    if (cin.fail()) {
        cin.clear();
        cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        cout << "Invalid input. Please enter an integer." << std::endl;
        mainMenu();
        return;
    }

    switch (k) {
        case 0:
            exit(0);
        case 1:
            printToyGraph();
            break;
        case 2:

            break;
        case 3:

            break;
        default:
            std::cout << "Option doesn't exist. Please insert a valid option." << std::endl;
            chooseGraph(); break;

    }
}

void Menu::printToyGraph() {
    std::cout << std::endl;
    std::cout << "#######################################################################" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##                 * * *   Toy Graphs   * * *                        ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      1 - Shipping Graph                                           ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      2 - Stadium Graph                                            ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      2 - Tourism Graph                                            ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      0-> Exit                                                     ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "#######################################################################" << std::endl << std::endl;


    int k;
    std::cout << "Option: ";
    std::cin >> k;

    if (cin.fail()) {
        cin.clear();
        cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        cout << "Invalid input. Please enter an integer." << std::endl;
        mainMenu();
        return;
    }

    switch (k) {
        case 0:
            exit(0);
        case 1:
            script.read_shipping();
            g = script.getShipGraph();
            mainMenu();
            break;
        case 2:
            script.read_stadiums();
            g = script.getStGraph();
            mainMenu();
            break;
        case 3:
            script.read_tourism();
            g = script.getTmGraph();
            mainMenu();
            break;
        default:
            std::cout << "Option doesn't exist. Please insert a valid option." << std::endl;
            printToyGraph(); break;

    }
}
