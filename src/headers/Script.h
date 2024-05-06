
#ifndef TRAVELLINGSALESMAN_SCRIPT_H
#define TRAVELLINGSALESMAN_SCRIPT_H

#include <istream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <cmath>

#include "Graph.h"

using namespace std;


class Script {
private:
    Graph<int> * shipGraph = new Graph<int>();
    Graph<int> * stGraph = new Graph<int>();
    Graph<int> * tmGraph = new Graph<int>();

    Graph<int> * rwg_g1 = new Graph<int>();
public:

    void read_shipping();
    void read_stadiums();
    void read_tourism();

    void read_rwg_g1();

    Graph<int> *getShipGraph() const;
    Graph<int> *getStGraph() const;
    Graph<int> *getTmGraph() const;

    Graph<int> *getRealWorldGraph1() const;
};


#endif //TRAVELLINGSALESMAN_SCRIPT_H
