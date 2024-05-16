#include "../headers/Script.h"

#include <iostream>

/**
 * @brief Reads shipping data from a CSV file and populates the shipping graph.
 * @details The file should contain edges in the format: origin, destination, distance.
 * Each edge is added in both directions to the graph.
 * @note The file path is hardcoded to "../datasets/toyGraphs/shipping.csv".
 */

void Script::read_shipping() {
    ifstream File1("../datasets/toyGraphs/shipping.csv");
    if (File1.is_open()){
        string line;
        double distance;
        int origin, destination;
        getline(File1, line);


        while(getline(File1, line)){
            istringstream iss(line);
            string orig, dest, dist;

            getline(iss, orig, ',');
            getline(iss, dest, ',');
            getline(iss, dist, ',');

            origin = stoi(orig);
            destination = stoi(dest);
            distance = stod(dist);

            // Add vertices to the graph with latitude and longitude parameters
            shipGraph->addVertex(origin, 0.0, 0.0);
            shipGraph->addVertex(destination, 0.0, 0.0);
            shipGraph->addEdge(origin, destination, distance);
            shipGraph->addEdge(destination, origin, distance);
        }
    }
    File1.close();
}

/**
 * @brief Reads stadium data from a CSV file and populates the stadium graph.
 * @details The file should contain edges in the format: origin, destination, distance.
 * Each edge is added in both directions to the graph.
 * @note The file path is hardcoded to "../datasets/toyGraphs/stadiums.csv".
 */

void Script::read_stadiums() {
    ifstream File1("../datasets/toyGraphs/stadiums.csv");
    if (File1.is_open()){
        string line;
        double distance;
        int origin, destination;
        getline(File1, line);

        while(getline(File1, line)){
            istringstream iss(line);
            string orig, dest, dist;

            getline(iss, orig, ',');
            getline(iss, dest, ',');
            getline(iss, dist, ',');

            origin = stoi(orig);
            destination = stoi(dest);
            distance = stod(dist);

            stGraph->addVertex(origin, 0.0, 0.0);
            stGraph->addVertex(destination, 0.0, 0.0);
            stGraph->addEdge(origin, destination, distance);
            stGraph->addEdge(destination, origin, distance);

        }
    }
    File1.close();
}

/**
 * @brief Reads tourism data from a CSV file and populates the tourism graph.
 * @details The file should contain edges in the format: origin, destination, distance, labelOrigin, labelDestination.
 * Each edge is added in both directions to the graph.
 * @note The file path is hardcoded to "../datasets/toyGraphs/tourism.csv".
 */

void Script::read_tourism() {
    ifstream File1("../datasets/toyGraphs/tourism.csv");
    if (File1.is_open()){
        string line;
        double distance;
        int origin, destination;
        getline(File1, line);

        while(getline(File1, line)){
            istringstream iss(line);
            string orig, dest, dist, labelOrig, labelDest;

            getline(iss, orig, ',');
            getline(iss, dest, ',');
            getline(iss, dist, ',');
            getline(iss, labelOrig, ',');
            getline(iss, labelDest, ',');

            origin = stoi(orig);
            destination = stoi(dest);
            distance = stod(dist);

            tmGraph->addVertex(origin, 0.0, 0.0);
            tmGraph->addVertex(destination, 0.0, 0.0);
            tmGraph->addEdge(origin, destination, distance);
            tmGraph->addEdge(destination, origin, distance);
        }
    }
    File1.close();
}

/**
 * @brief Reads real-world graph 1 data from CSV files and populates the graph.
 * @details Reads node coordinates from "nodes.csv" and edges from "edges.csv".
 * Each edge is added in both directions to the graph.
 * @note The file paths are hardcoded to "../datasets/Real-World Graphs/graph1/nodes.csv" and "../datasets/Real-World Graphs/graph1/edges.csv".
 */

void Script::read_rwg_g1() {

    ifstream File2("../datasets/Real-World Graphs/graph1/nodes.csv");
    unordered_map<int, pair<double,double>> coordinates;

    if (File2.is_open()) {
        string line2;
        int id;
        double longitude, latitude;
        getline(File2, line2);

        while (getline(File2, line2)) {
            istringstream iss2(line2);
            string i, l, ll;

            getline(iss2, i, ',');
            getline(iss2, l, ',');
            getline(iss2, ll, ',');

            id = stoi(i);
            longitude = stod(l);
            latitude = stod(ll);

            coordinates[id] = {longitude, latitude};
        }
    }
    File2.close();

    ifstream File1("../datasets/Real-World Graphs/graph1/edges.csv");

    if (File1.is_open()) {
        string line1;
        int origin, destination;
        double distance;
        getline(File1, line1);

        while (getline(File1, line1)) {
            istringstream iss1(line1);
            string orig, dest, dist;

            getline(iss1, orig, ',');
            getline(iss1, dest, ',');
            getline(iss1, dist, ',');

            origin = stoi(orig);
            destination = stoi(dest);
            distance = stod(dist);

            rwg_g1->addVertex(origin, coordinates[origin].first, coordinates[origin].second);
            rwg_g1->addVertex(destination, coordinates[destination].first, coordinates[destination].second);
            rwg_g1->addEdge(origin, destination, distance);
            rwg_g1->addEdge(destination, origin, distance);
        }

    }
    File1.close();

}

/**
 * @brief Reads extra fully connected graph data from a CSV file and populates the graph.
 * @details The file should contain edges in the format: origin, destination, distance.
 * Each edge is added in both directions to the graph.
 * @note The file path is hardcoded to "../datasets/Extra_Fully_Connected_Graphs/edges_25.csv".
 */

void Script::read_efcg_25() {
    ifstream File1("../datasets/Extra_Fully_Connected_Graphs/edges_25.csv");

    if (File1.is_open()) {
        string line1;
        int origin, destination;
        double distance;
        getline(File1, line1);

        while (getline(File1, line1)) {
            istringstream iss1(line1);
            string orig, dest, dist;

            getline(iss1, orig, ',');
            getline(iss1, dest, ',');
            getline(iss1, dist, ',');

            origin = stoi(orig);
            destination = stoi(dest);
            distance = stod(dist);

            efcg_25->addVertex(origin, 0.0, 0.0);
            efcg_25->addVertex(destination, 0.0, 0.0);
            efcg_25->addEdge(origin, destination, distance);
            efcg_25->addEdge(destination, origin, distance);
        }

    }
    File1.close();
}

/**
 * @brief Gets the shipping graph.
 * @return A pointer to the shipping graph.
 */

Graph<int> * Script::getShipGraph() const {
    return shipGraph;
}

/**
 * @brief Gets the stadium graph.
 * @return A pointer to the stadium graph.
 */


Graph<int> * Script::getStGraph() const {
    return stGraph;
}

/**
 * @brief Gets the tourism graph.
 * @return A pointer to the tourism graph.
 */

Graph<int> * Script::getTmGraph() const {
    return tmGraph;
}

/**
 * @brief Gets the real-world graph 1.
 * @return A pointer to the real-world graph 1.
 */

Graph<int> *Script::getRealWorldGraph1() const {
    return rwg_g1;
}

/**
 * @brief Gets the extra fully connected graph with 25 nodes.
 * @return A pointer to the extra fully connected graph.
 */

Graph<int> *Script::getExtraFulllyConnected25() const {
    return efcg_25;
}
