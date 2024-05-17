

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
            tspSolver.uniteAllClusterTours(g, tspSolver.kMeansClustering(g));
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
            printExtraFullyConnected();
            break;
        case 3:
            printRealWorldGraph();
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

/**
 * @brief Displays the real world graph selection menu and processes user input to select a specific real world graph.
 * @details This function shows a menu with different toy graphs for the user to choose from.
 * Depending on the user input, it loads the selected real world graph by calling the appropriate function from the Script class.
 */

void Menu::printRealWorldGraph() {
    std::cout << std::endl;
    std::cout << "#######################################################################" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##                * * *   Real World Graphs   * * *                  ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      1 - Graph 1                                                  ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      2 - Graph 2                                                  ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      3 - Graph 3                                                  ##" << std::endl;
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
            script.read_rwg("../cmake-build-debug/datasets/Real-World Graphs/graph1/nodes.csv", "../cmake-build-debug/datasets/Real-World Graphs/graph1/edges.csv");
            g = script.getRealWorldGraph();
            mainMenu();
            break;
        case 2:
            script.read_rwg("../cmake-build-debug/datasets/Real-World Graphs/graph2/nodes.csv", "../cmake-build-debug/datasets/Real-World Graphs/graph2/edges.csv");
            g = script.getRealWorldGraph();
            mainMenu();
            break;
        case 3:
            script.read_rwg("../cmake-build-debug/datasets/Real-World Graphs/graph2/nodes.csv", "../cmake-build-debug/datasets/Real-World Graphs/graph2/edges.csv");
            g = script.getRealWorldGraph();
            mainMenu();
            break;
        default:
            std::cout << "Option doesn't exist. Please insert a valid option." << std::endl;
            printRealWorldGraph();
            break;
    }
}

/**
 * @brief Displays the Extra Fully Connected selection menu and processes user input to select a specific Extra Fully Connected graph.
 * @details This function shows a menu with different tExtra Fully Connected graphs for the user to choose from.
 * Depending on the user input, it loads the selected Extra Fully Connected graph by calling the appropriate function from the Script class.
 */

void Menu::printExtraFullyConnected() {
    std::cout << std::endl;
    std::cout << "#######################################################################" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##         * * *   Extra Fully Connected  Graphs  * * *              ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      1 - 25 Edges                                                 ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      2 - 50 Edges                                                 ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      3 - 75 Edges                                                 ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      4 - 100 Edges                                                ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      5 - 200 Edges                                                ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      6 - 300 Edges                                                ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      7 - 400 Edges                                                ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      8 - 500 Edges                                                ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      9 - 600 Edges                                                ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      10 - 700 Edges                                               ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      11 - 800 Edges                                               ##" << std::endl;
    std::cout << "##                                                                   ##" << std::endl;
    std::cout << "##      12 - 900 Edges                                               ##" << std::endl;
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
            script.read_efcg("../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/nodes.csv", "../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/edges_25.csv");
            g = script.getExtraFullyConnected();
            mainMenu();
            break;
        case 2:
            script.read_efcg("../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/nodes.csv", "../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/edges_50.csv");
            g = script.getExtraFullyConnected();
            mainMenu();
            break;
        case 3:
            script.read_efcg("../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/nodes.csv", "../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/edges_75.csv");
            g = script.getExtraFullyConnected();
            mainMenu();
            break;
        case 4:
            script.read_efcg("../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/nodes.csv", "../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/edges_100.csv");
            g = script.getExtraFullyConnected();
            mainMenu();
            break;
        case 5:
            script.read_efcg("../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/nodes.csv", "../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/edges_200.csv");
            g = script.getExtraFullyConnected();
            mainMenu();
            break;
        case 6:
            script.read_efcg("../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/nodes.csv", "../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/edges_300.csv");
            g = script.getExtraFullyConnected();
            mainMenu();
            break;
        case 7:
            script.read_efcg("../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/nodes.csv", "../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/edges_400.csv");
            g = script.getExtraFullyConnected();
            mainMenu();
            break;
        case 8:
            script.read_efcg("../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/nodes.csv", "../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/edges_500.csv");
            g = script.getExtraFullyConnected();
            mainMenu();
            break;
        case 9:
            script.read_efcg("../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/nodes.csv", "../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/edges_600.csv");
            g = script.getExtraFullyConnected();
            mainMenu();
            break;
        case 10:
            script.read_efcg("../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/nodes.csv", "../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/edges_700.csv");
            g = script.getExtraFullyConnected();
            mainMenu();
            break;
        case 11:
            script.read_efcg("../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/nodes.csv", "../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/edges_800.csv");
            g = script.getExtraFullyConnected();
            mainMenu();
            break;
        case 12:
            script.read_efcg("../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/nodes.csv", "../cmake-build-debug/datasets/Extra_Fully_Connected_Graphs/edges_900.csv");
            g = script.getExtraFullyConnected();
            mainMenu();
            break;
        default:
            std::cout << "Option doesn't exist. Please insert a valid option." << std::endl;
            printExtraFullyConnected();
            break;
    }
}
