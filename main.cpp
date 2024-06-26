#include "capd/capdlib.h"

#include "main.h"
#include "solvers.h"
#include "neuron.h"
#include "capd_error.h"
#include "print_state.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using namespace std;
using namespace capd;



Statistics statistics;

void update_statistics(int bound, int max_number)
{
    if (bound == NO_SOLUTION )
        statistics.no_solution++;
    else if ( bound == ERROR_MAX_DERIVATIVE )
        statistics.max_derivative++;
    else if ( bound == ERROR_MAX_SUBDIVISIONS )
        statistics.max_subdivisions++;
    else if ( bound > max_number )
        statistics.max_number++;
}

int iterate_tolerance(IMap& target, IVector x, IVector p, int max_derivative, int max_number, double tolerance)
{
    const double init = 1.1e-1;
    assert(tolerance < init);
    int bound;
    for (double delta = init; delta > tolerance; delta /= 2.) {
        bound = bifurcation_order(target, x, p, max_derivative, delta);
        if(bound < 0) {
            continue;
        }
        if(bound > max_number) {
            continue;
        }
        return(bound);
    }
    return(bound);
}

void bifurcation_order_wrapper(IMap& target, State& state, IVector x, int max_number, int max_derivative)
{
    /**
        Sorts through intervals in state.special. Either moves them to state.regular, state.verified, or leaves them in state.special, depending on the output of bifurcation_order. Does not subdivide.
    **/

    vector<IVector> new_special;
    int bound;

    for_each (state.special.begin(), state.special.end(), [&](IVector p) {
        bound = iterate_tolerance(target, x, p, max_derivative, max_number, state.tolerance);

        update_statistics(bound, max_number);

        if (bound == NO_SOLUTION ) {
            state.regular.push_back(p);
        } else if ( bound == ERROR_MAX_DERIVATIVE || bound == ERROR_MAX_SUBDIVISIONS || bound > max_number ) {
            new_special.push_back(p);
        } else {
            state.verified.push_back(p);
        }

        cout << '\r';
        cout << setw(10) << right << state.regular.size()
             << setw(10) << right << state.special.size()
             << setw(10) << right << state.verified.size()
             << setw(13) << right << new_special.size() << "  "
             << setw(15) << left << state.tolerance << flush;
    });
    cout << endl;

    state.special.clear();
    state.special = new_special;
}

int newton_wrapper(IMap& target, vector<vector<IVector>>& regular_components, int max_number, IVector x, bool verbose)
{
    /**
        Verifies the conjectured maximum number of solutions, max_number, for a single point in each component in regular_components.
    **/

    int result;
    IVector x_newton, p;

    for (auto & comp : regular_components)
    {
        x_newton = x;

        p = comp[0];
        p = midVector(p);

        while ( (result = newton_method(target, x_newton, p)) == NEWTON_DOUBLE_COUNT )
        {
            if(verbose) cout << "Double counted solution" << endl;
            x_newton[0] += interval(0, 1e-5);
            if(verbose) cout << "\\x = " << x_newton << " and \\p = " << p << endl;
        }

        if ( result == NEWTON_EPSILON_WIDTH ) {
            if(verbose) {
                cout << "Couldn't resolve solutions for parameters " << p
                << " - try a smaller inflation or different parameters" << endl;
            }
            return(NEWTON_EPSILON_WIDTH);
        } else if ( result > max_number ) {
            if(verbose) {
                cout << "Newton-verified number of solutions for parameters " << p << " is " << result << ". Terminating." << endl;
            }
            return(-3);
        }
    }
    return(0);
}

void interval_hull(vector<IVector>& intervals, IVector* hull)
{
    assert(!intervals.empty());
    IVector current_hull = intervals[0];
    for(auto & vec : intervals) {
        current_hull = intervalHull(current_hull, vec);
    }
    *hull = current_hull;
}

int prove_bound(IMap& target, IVector x, IVector p, int max_number, int max_derivative, double tolerance)
{
    State state;
    state.special.push_back(p);
    state.tolerance = tolerance;

    cout << setw(10) << right << "regular"
         << setw(10) << right << "special"
         << setw(10) << right << "verified"
         << setw(13) << right << "new_special" << "  "
         << setw(15) << left << "tolerance" << endl;

    bool print_params = false;

    do {
        bisection(target, x, state); // subdivides special boxes
        bifurcation_order_wrapper(target, state, x, max_number, max_derivative);
        if(state.tolerance < 1e-4) {print_params = true;}////////////////////
        if(print_params) {
            IVector hull;
            interval_hull(state.special, &hull);
            std::cout << "Interval hull of special is " << hull << std::endl;
        }////////////////////
        if(print_params) {
            cout << "DIAGNOSTIC!! IGNORE SUCCESS!!" << endl;
            break;
        }////////////////////
        state.tolerance /= 2.;
    } while (!state.special.empty() && state.tolerance > PROVE_TOLERANCE);

    if(p.dimension() == 2) {
        print_state(state);
    }

    if (!state.special.empty())
        return(ERROR_PROVE_TOLERANCE);

    vector<vector<IVector>> regular_components;
    find_connected_components(state.regular, regular_components);
    cout << "\nthere are " << regular_components.size() << " regular components" << endl;
    return(newton_wrapper(target, regular_components, max_number, x, false));
}

int main(int argc, char **argv)
{
    int max_derivative = 3; // maximum multiplicity
    int max_number = 3;

    int phase_dim = 1;
    int param_dim = 2;

    IMap target = get_target(max_derivative, phase_dim, param_dim);

    IVector x(1), p(2);

    x[0] = interval(-2, 2 + 1e-1);//(0.2101, 0.22);

    vector<double> bounds(argc - 1);

    for (int i = 0; i < argc - 1; i++) {
        bounds[i] = std::strtod(argv[i + 1], nullptr);
    }

    p[0] = interval(bounds[0], bounds[1]);///(0.65, 0.7);
    p[1] = interval(bounds[2], bounds[3]);///(0.65, 0.7);

    double tolerance = 1e-1;

    auto t1 = std::chrono::high_resolution_clock::now();
    int result = prove_bound(target, x, p, max_number, max_derivative, tolerance);
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "\nprove_bound took "
        << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
        << " microseconds" << endl;

    if (result == 0) {
        cout << "\n\nSuccess! Maximum of " << max_number << " fixed points in " << p << "." << endl;
    } else {
        cout << "Failed with error code " << result << "." << endl;
    }

    cout << "\nStatistics:" << endl;
    cout << "no_solution = " << statistics.no_solution << endl;
    cout << "max_derivative = " << statistics.max_derivative << endl;
    cout << "max_subdivisions = " << statistics.max_subdivisions << endl;
    cout << "max_number = " << statistics.max_number << endl;
    cout << "components_dur in microseconds = " << statistics.components_dur.count() << endl;
    cout << "refine_dur in microseconds = " << statistics.refine_dur.count() << endl;
    cout << "order_dur in microseconds = " << statistics.order_dur.count() << endl;

    return(0);
}
