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



Statistics statistics;

void bifurcation_order_wrapper(IMap& target, State& state, int max_number, int max_derivative)
{

    /**
        Sorts through intervals in state.special. Either moves them to state.regular, state.verified, or leaves them in state.special, depending on the output of bifurcation_order. Does not subdivide.
    **/

    vector<SpecialBox> new_special;

    int result;

    for_each (state.special.begin(), state.special.end(), [&](SpecialBox s_box) {

        cout << '\r';
        cout << setw(10) << right << state.regular.size()
             << setw(10) << right << state.special.size()
             << setw(10) << right << state.verified.size()
             << setw(13) << right << new_special.size() << "  "
             << setw(15) << left << state.tolerance << flush;

        result = bifurcation_order(target, s_box, max_derivative);

        if (result == ERROR_MAX_DERIVATIVE) {
            statistics.max_derivative++;
        } else if ( result > max_number ) {
            statistics.max_number++;
        }

        ///if ( result == ERROR_NO_BIFURCATION )

        ///    state.regular.push_back(s_box.parameters);

        if ( result == ERROR_MAX_DERIVATIVE || result > max_number )

            new_special.push_back(s_box);

        else

            state.verified.push_back(s_box.parameters);

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

int automatic(IMap& target, const IVector& phase_space, IVector& parameter_space, int max_number, int max_derivative, double tolerance)
{

    State state;

    SpecialBox s_box;
    s_box.parameters = parameter_space;
    s_box.phase_boxes.push_back(phase_space);
    state.special.push_back(s_box);

    state.tolerance = tolerance;

    cout << setw(10) << right << "regular"
         << setw(10) << right << "special"
         << setw(10) << right << "verified"
         << setw(13) << right << "new_special" << "  "
         << setw(15) << left << "tolerance" << endl;

    while (!state.special.empty() && state.tolerance > AUTO_TOLERANCE) {
        bisection(target, state); // subdivides special boxes
        bifurcation_order_wrapper(target, state, max_number, max_derivative);
        state.tolerance /= 2.;
    }

    if (!state.special.empty())
        return(ERROR_AUTO_TOLERANCE);///should really write special boxes to disk

    vector<vector<IVector>> regular_components;
    find_connected_components(state.regular, regular_components);
    return(newton_wrapper(target, regular_components, max_number, phase_space, false));

}

int main()
{
    int max_derivative = 3; // maximum multiplicity
    int max_number = 5;

    IMap target = get_target(max_derivative);

    IVector x(1), p(3);

    x[0] = interval(-1.0, 1.0);///interval(0.2101, 0.22);

    interval diameter = interval(-1,1) * 1e-2;
    p[0] = interval(-1.) + diameter;
    p[1] = interval(2.) + diameter;
    p[2] = interval(1.) + diameter;

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