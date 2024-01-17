using namespace capd;
using namespace std;

void bisection(IMap& target, IVector x, IVector p, vector<IVector>& regular, vector<IVector>& special, double tolerance);

int newton_method(IMap& target, IVector x, IVector p);

int bifurcation_order(IMap& target, IVector x, IVector p, int max_derivative, double tolerance);
