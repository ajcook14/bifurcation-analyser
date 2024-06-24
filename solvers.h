#ifndef SOLVERS_H
#define SOLVERS_H
#include <chrono>

using namespace capd;
using namespace std;

struct State {
    vector<IVector> regular;    // contains only intervals in parameter space with no bifurcations
    vector<IVector> special;    // contains intervals in parameter space that might contain bifurcations
    vector<IVector> verified;   // intervals verified to have no more than max_number fixed points
    double tolerance;
};

struct Statistics {
    int no_solution = 0;
    int max_derivative = 0;
    int max_subdivisions = 0;
    int max_number = 0;

    std::chrono::duration<long, std::ratio<1, 1000000>> components_dur = std::chrono::microseconds(0);
    std::chrono::duration<long, std::ratio<1, 1000000>> refine_dur = std::chrono::microseconds(0);
    std::chrono::duration<long, std::ratio<1, 1000000>> order_dur = std::chrono::microseconds(0);
};

extern bool print_params;///////////////

void bisection(IMap& target, IVector x, State& state);
int newton_method(IMap& target, IVector x, IVector p);
void find_connected_components(vector<IVector>& intervals, vector<vector<IVector>>& components);
int bifurcation_order(IMap& target, IVector x, IVector p, int max_derivative, double tolerance);

#endif // SOLVERS_H
