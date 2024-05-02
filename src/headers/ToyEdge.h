
#ifndef TRAVELLINGSALESMAN_TOYEDGE_H
#define TRAVELLINGSALESMAN_TOYEDGE_H

#include <string>
using namespace std;
class ToyEdge{
private:
    string orig;
    string dest;
    string type;
    double distance;

public:
    ToyEdge(string orig, string dest, string type, double distance);
    const string getOrig() const;
    const int getId() const;
    const string getDest() const;
    const double getDistance() const;
    const string getType() const;
};

#endif //TRAVELLINGSALESMAN_STADIUM_H
