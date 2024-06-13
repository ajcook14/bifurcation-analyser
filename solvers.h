
#ifndef BIFURCATION_ANALYSER_SOLVERS_H
#define BIFURCATION_ANALYSER_SOLVERS_H

#endif //BIFURCATION_ANALYSER_SOLVERS_H

using namespace capd;
using namespace std;

struct SpecialBox {

    IVector parameters;             // box of parameters

    vector<IVector> phase_boxes;    // collection of boxes containing all solutions in phase space (if any)

};

struct State {

    vector<IVector> regular;    // contains only intervals in parameter space with no bifurcations

    vector<SpecialBox> special; // all intervals in parameter space that might contain bifurcations

    vector<IVector> verified;   // intervals verified to have no more than max_number fixed points

    double tolerance;

};

struct Statistics {

    int no_bifurcation = 0;
    int max_derivative = 0;
    int max_subdivisions = 0;
    int max_number = 0;

};

void bisection(IMap&, State&);

int newton_method(IMap&, IVector x, IVector p);

void find_connected_components(vector<IVector>&, vector<vector<IVector>>&);

int bifurcation_order(IMap&, SpecialBox, int);

int refine_measure(IMap&/******/, SpecialBox*, double);