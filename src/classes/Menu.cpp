

#include "../headers/Menu.h"

/**
 * @brief Constructs a Menu object with a TSPSolver and a Script.
 * @param tspSolver An instance of the TSPSolver class.
 * @param script An instance of the Script class.
 */
Menu::Menu(TSPSolver tspSolver, Script script) : tspSolver(tspSolver), script(script) {}

/**
 * @brief Displays the main menu and processes user input to select TSP algorithms or other options.
 * @details This function shows a menu with several TSP algorithms and options for the user to choose from.
 * Depending on the user input, it calls the appropriate function from the TSPSolver or navigates to a sub-menu.
 */
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
            tspSolver.calculateTSP(g);
            mainMenu();
            break;
        case 2:
            tspSolver.calculateTriangleTSP(g);
            mainMenu();
            break;
        case 3:
            mainMenu();
            break;
        case 4:
            tspSolver.calculateNearestNeighborTSP(g);
            mainMenu();
            break;
        case 5:
            chooseGraph();
            break;
        default:
            std::cout << "Please insert a valid option." << std::endl;
            mainMenu();
            break;
    }
}

/**
 * @brief Displays the graph selection menu and processes user input to select the type of graph.
 * @details This function shows a menu with different types of graphs for the user to choose from.
 * Depending on the user input, it calls the appropriate function to load the selected graph.
 */
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
    std::cout << "##      3 - Real World Graphs                                        ##" << std::endl;
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
            // Code to handle Extra Fully Connected Graphs selection
            break;
        case 3:
            // Code to handle Real World Graphs selection
            break;
        default:
            std::cout << "Option doesn't exist. Please insert a valid option." << std::endl;
            chooseGraph();
            break;
    }
}

/**
 * @brief Displays the toy graph selection menu and processes user input to select a specific toy graph.
 * @details This function shows a menu with different toy graphs for the user to choose from.
 * Depending on the user input, it loads the selected toy graph by calling the appropriate function from the Script class.
 */
 
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
    std::cout << "##      3 - Tourism Graph                                            ##" << std::endl;
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
            printToyGraph();
            break;
    }
}
