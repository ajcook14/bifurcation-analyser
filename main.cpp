#include "capd/capdlib.h"

#include "main.h"
#include "solvers.h"
#include "neuron.h"
#include "capd_error.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace capd;



void update_statistics(int bound, int max_number)
{

    if ( bound == ERROR_NO_BIFURCATION )
        statistics.no_bifurcation++;
    else if ( bound == ERROR_MAX_DERIVATIVE )
        statistics.max_derivative++;
    else if ( bound == ERROR_MAX_SUBDIVISIONS )
        statistics.max_subdivisions++;
    else if ( bound > max_number )
        statistics.max_number++;

}

void bifurcation_order_wrapper(IMap& target, State& state, IVector x, int max_number, int max_derivative)
{

    /**
        Sorts through intervals in state.special. Either moves them to state.regular, state.verified, or leaves them in state.special, depending on the output of bifurcation_order. Does not subdivide.
    **/

    vector<IVector> new_special;

    int bound;

    for_each (state.special.begin(), state.special.end(), [&](IVector p) {

        cout << '\r';
        cout << setw(10) << right << state.regular.size()
        << setw(10) << right << state.special.size()
        << setw(10) << right << state.verified.size()
        << setw(10) << right << new_special.size()
        << setw(15) << right << state.tolerance << flush;

        bound = bifurcation_order(target, x, p, max_derivative, state.tolerance);

        update_statistics(bound, max_number);

        if ( bound == ERROR_NO_BIFURCATION )

            state.regular.push_back(p);

        else if ( bound == ERROR_MAX_DERIVATIVE || bound == ERROR_MAX_SUBDIVISIONS || bound > max_number )

            new_special.push_back(p);

        else

            state.verified.push_back(p);

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

    for (int i = 0; i < regular_components.size(); i++)
    {

        x_newton = x;

        p = regular_components[i][0];
        p = midVector(p);

        while ( (result = newton_method(target, x_newton, p)) == NEWTON_DOUBLE_COUNT )
        {

            if(verbose) cout << "Double counted solution" << endl;

            x_newton[0] += interval(0, 1e-5);

            if(verbose) cout << "\\x = " << x_newton << " and \\p = " << p << endl;

        }

        if ( result == NEWTON_EPSILON_WIDTH ) {

            if(verbose) cout << "Couldn't resolve solutions - try a smaller inflation or different parameters" << endl;

            return(NEWTON_EPSILON_WIDTH);

        } else if ( result > max_number ) {

            if(verbose) cout << "Newton-verified number of solutions for parameters " << p << " is " << result << ". Terminating." << endl;

            return(-3);

        }

    }

    return(0);

}

int automatic(IMap& target, IVector x, IVector p, int max_number, int max_derivative, double tolerance)
{

    State state;

    state.special.push_back(p);

    state.tolerance = tolerance;

    cout << setw(10) << right << "regular"
    << setw(10) << right << "special"
    << setw(10) << right << "verified"
    << setw(10) << right << "new_special"
    << setw(15) << right << "tolerance" << endl;

    do {

        bisection(target, x, state); // subdivides special boxes

        bifurcation_order_wrapper(target, state, x, max_number, max_derivative);


        state.tolerance /= 2.;

    } while (!state.special.empty() && state.tolerance > AUTO_TOLERANCE);

    if (!state.special.empty())

        return(ERROR_AUTO_TOLERANCE);

    vector<vector<IVector>> regular_components;

    find_connected_components(state.regular, regular_components);

    return(newton_wrapper(target, regular_components, max_number, x, false));

}

Statistics statistics;

int main()
{

    int max_derivative = 3; // maximum multiplicity

    int max_number = 5;

    IMap target = get_target(max_derivative);

    IVector x(1), p(3);

    x[0] = interval(-1, 1);//(0.2101, 0.22);

    p[0] = interval(-1.) + interval(-0.0, 0.5);
    p[1] = interval(2.) + interval(-0.5, 0.5);
    p[2] = interval(1.) + interval(-0.5, 0.5);

    double tolerance = 1e-1;

    int result = automatic(target, x, p, max_number, max_derivative, tolerance);

    if (result == 0)
        cout << "\n\nSuccess! Maximum of " << max_number << " fixed points in " << p << "." << endl;
    else
        cout << "Failed with error code " << result << "." << endl;

    cout << "\nStatistics:" << endl;
    cout << "no_bifurcation = " << statistics.no_bifurcation << endl;
    cout << "max_derivative = " << statistics.max_derivative << endl;
    cout << "max_subdivisions = " << statistics.max_subdivisions << endl;
    cout << "max_number = " << statistics.max_number << endl;

    return(0);

}
