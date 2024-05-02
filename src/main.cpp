#include <iostream>

#include "headers/Script.h"

using namespace std;

int main() {
    Script script;
    script.read_tourism();
    for (Vertex<int>* v: script.getTmGraph()->getVertexSet()){
        for (Edge<int> * e: v->getAdj()) {
            cout << e->getOrig()->getInfo() << " " << e->getDest()->getInfo() << " " <<  e->getWeight()  << endl;
        }
    }
    return 0;
}
