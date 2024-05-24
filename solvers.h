using namespace capd;
using namespace std;

struct State {

    vector<IVector> regular;    // contains only intervals in parameter space with no bifurcations

    vector<IVector> special;    // contains intervals in parameter space that might contain bifurcations

    vector<IVector> verified;   // intervals verified to have no more than max_number fixed points

    double tolerance;

};

void bisection(IMap& target, IVector x, State& state);

int newton_method(IMap& target, IVector x, IVector p);

void find_connected_components(vector<IVector>& intervals, vector<vector<IVector>>& components);

int bifurcation_order(IMap& target, IVector x, IVector p, int max_derivative, double tolerance);

