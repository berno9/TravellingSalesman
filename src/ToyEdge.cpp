
#include "ToyEdge.h"

ToyEdge::ToyEdge(string orig, string dest, string type, double distance) :
        orig(orig), dest(dest),  type(type), distance(distance) {}

const string ToyEdge::getOrig() const { return orig; }

const string ToyEdge::getDest() const { return dest; }

const double ToyEdge::getDistance() const { return distance; }

const string ToyEdge::getType() const { return type; }

