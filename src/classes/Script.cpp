#include "../headers/Script.h"

#include <iostream>

/**
 * @brief Reads shipping data from a CSV file and populates the shipping graph.
 * @details The file should contain edges in the format: origin, destination, distance.
 * Each edge is added in both directions to the graph.
 * @note The file path is hardcoded to "../datasets/toyGraphs/shipping.csv".
 */

void Script::read_shipping() {
    ifstream File1("../cmake-build-debug/datasets/toyGraphs/shipping.csv");
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
    ifstream File1("../cmake-build-debug/datasets/toyGraphs/stadiums.csv");
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
    ifstream File1("../cmake-build-debug/datasets/toyGraphs/tourism.csv");
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
 * @brief Reads real-world graph  data from CSV files and populates the graph.
 * @details Reads node coordinates from "nodes.csv" and edges from "edges.csv".
 * Each edge is added in both directions to the graph.
 */

void Script::read_rwg(string s1, string s2) {

    ifstream File2(s1);
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

    ifstream File1(s2);

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

            rwg->addVertex(origin, coordinates[origin].first, coordinates[origin].second);
            rwg->addVertex(destination, coordinates[destination].first, coordinates[destination].second);
            rwg->addEdge(origin, destination, distance);
            rwg->addEdge(destination, origin, distance);
        }

    }
    File1.close();

}

/**
 * @brief Reads extra fully connected graph data from a CSV file and populates the graph.
 * @details The file should contain edges in the format: origin, destination, distance.
 * Each edge is added in both directions to the graph.
 */

void Script::read_efcg(string s1, string s2) {
    ifstream File2(s1);
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

    ifstream File1(s2);

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

            efcg->addVertex(origin, coordinates[origin].first, coordinates[origin].second);
            efcg->addVertex(destination, coordinates[destination].first, coordinates[destination].second);
            efcg->addEdge(origin, destination, distance);
            efcg->addEdge(destination, origin, distance);
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
 * @brief Gets a real-world graph.
 * @return A pointer to a real-world graph.
 */

Graph<int> *Script::getRealWorldGraph() const {
    return rwg;
}

/**
 * @brief Gets the extra fully connected graph with a required number of nodes.
 * @return A pointer to the extra fully connected graph.
 */

Graph<int> *Script::getExtraFullyConnected() const {
    return efcg;
}
