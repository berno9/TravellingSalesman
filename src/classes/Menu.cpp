
#include "../headers/Menu.h"

Menu::Menu(TSPSolver tspSolver, Script script) : tspSolver(tspSolver), script(script) {}

void Menu::mainMenu() {
    std::cout << std::endl;
    std::cout << "#######################################################################" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##                    * * *   Consult   * * *                        ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      1 - TSP with Backtracking                                    ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      2 - TSP with Triangular Approximation Heuristic              ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
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
            
            break;
        case 2:

        default:
            std::cout << "Option doesn't exist. Please insert a valid option." << std::endl;
            mainMenu(); break;

    }
}

void Menu::switchGraph() {
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
            mainMenu(); break;

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

            break;
        case 2:

            break;
        case 3:

            break;
        default:
            std::cout << "Option doesn't exist. Please insert a valid option." << std::endl;
            mainMenu(); break;

    }
}

void Menu::solutionTSPBack(Graph<int>* g) {

}