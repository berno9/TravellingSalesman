
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

            // Add vertices to the graph with latitude and longitude parameters
            shipGraph->addVertex(origin, 0.0, 0.0);
            shipGraph->addVertex(destination, 0.0, 0.0);
            shipGraph->addEdge(origin, destination, distance);
            bool dir = true;
            for (Vertex<int>* v: shipGraph->getVertexSet()){
                if (v->getInfo() == destination){
                    for (Edge<int>* e: v->getAdj()){
                        if (e->getDest()->getInfo() == origin) {
                            dir = false;
                        }
                    }
                }
            }

            if (dir) {
                shipGraph->addEdge(destination, origin, distance);
            }
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

            stGraph->addVertex(origin, 0.0, 0.0);
            stGraph->addVertex(destination, 0.0, 0.0);
            stGraph->addEdge(origin, destination, distance);

            bool dir = true;
            for (Vertex<int>* v: stGraph->getVertexSet()){
                if (v->getInfo() == destination){
                    for (Edge<int>* e: v->getAdj()){
                        if (e->getDest()->getInfo() == origin) {
                            dir = false;
                        }
                    }
                }
            }

            if (dir) {
                stGraph->addEdge(destination, origin, distance);
            }
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

            tmGraph->addVertex(origin, 0.0, 0.0);
            tmGraph->addVertex(destination, 0.0, 0.0);
            tmGraph->addEdge(origin, destination, distance);

            bool dir = true;
            for (Vertex<int>* v: tmGraph->getVertexSet()){
                if (v->getInfo() == destination){
                    for (Edge<int>* e: v->getAdj()){
                        if (e->getDest()->getInfo() == origin) {
                            dir = false;
                        }
                    }
                }
            }

            if (dir) {
                tmGraph->addEdge(destination, origin, distance);
            }
        }
    }
    File1.close();
}

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

Graph<int> *Script::getRealWorldGraph1() const {
    return rwg_g1;
}