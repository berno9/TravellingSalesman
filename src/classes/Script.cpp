
#include "../headers/Script.h"

#include <iostream>

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

            shipGraph->addVertex(origin);
            shipGraph->addVertex(destination);
            shipGraph->addEdge(origin, destination, distance);
        }
    }
    File1.close();
}

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

            stGraph->addVertex(origin);
            stGraph->addVertex(destination);
            stGraph->addEdge(origin, destination, distance);
        }
    }
    File1.close();
}

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

            tmGraph->addVertex(origin);
            tmGraph->addVertex(destination);
            tmGraph->addEdge(origin, destination, distance);
        }
    }
    File1.close();
}

Graph<int> * Script::getShipGraph() const {
    return shipGraph;
}

Graph<int> * Script::getStGraph() const {
    return stGraph;
}

Graph<int> * Script::getTmGraph() const {
    return tmGraph;
}