#include "capd/capdlib.h"

#include "solvers.h"
#include "neuron.h"
#include "capd_error.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace capd;



void write_vector(IVector& iv, ofstream* to)
// iv is dimension 2
{

    *to << iv[0].leftBound() << ',' << iv[0].rightBound() << ','
       << iv[1].leftBound() << ',' << iv[1].rightBound() << endl;

}

void find_connected_components(vector<IVector>& intervals, vector<vector<IVector>>& components)
{

    /**
        Writes the connected components of intervals to components.
        The vector intervals will get "sorted" into connected components.
    **/

    int n = intervals.size();

    int a[n];

    int c_index = 1;

    int k = 0;  // end index of current component

    for (int i = 0; i < n; i++) {

        a[i] = c_index;

        for (int j = k + 1; j < n; j++) {

            if ( !intersectionIsEmpty( intervals[i], intervals[j] ) ) {

                k++;

                std::swap( intervals[k], intervals[j] );

            }

        }

        if (k == i) {

            if (i < n) {

                c_index++;

                k++;

            }

            else break;

        }

    }

    // now write to components

    c_index = 1;

    vector<IVector> component;

    for (int i = 0; i < n; i++) {

        if (a[i] > c_index) {

            components.push_back(component);

            component.clear();

            c_index++;

        }

        component.push_back(intervals[i]);

    }

    if ( n > 0 ) components.push_back(component);

}

void bifurcation_order_wrapper(IMap& target, State& state, IVector x, int max_number, int max_derivative)
{

    /**
        Sorts through intervals in state.special. Either moves them to state.regular, state.verified, or leaves them in state.special, depending on the output of bifurcation_order. Does not subdivide.
    **/

    vector<IVector> new_special;

    int bound;

    for_each (state.special.begin(), state.special.end(), [&](IVector p) {

        bound = bifurcation_order(target, x, p, max_derivative, state.tolerance);

        if ( bound == ERROR_NO_BIFURCATION )

            state.regular.push_back(p);

        else if ( bound == ERROR_MAX_DERIVATIVE || bound > max_number )

            new_special.push_back(p);

        else

            state.verified.push_back(p);

    });

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

/*
void manual(IMap& target, int max_number, int max_derivative, IVector x, IVector p, double tolerance)
{

    vector<IVector>::iterator intvec_it;



    cout << "\nBISECTION:" << endl;

    vector<IVector> regular, special;

    cout << "\\x = " << x << " and \\p =  " << p << endl;
    bisection(target, x, p, regular, special, tolerance);

    cout << "Number of eliminated boxes = " << regular.size() << endl;
    cout << "Number of special boxes = " << special.size() << endl;

    vector<vector<IVector>> regular_components;

    find_connected_components(regular, regular_components);

    cout << "Number of components is " << regular_components.size() << endl;



    cout << "\nNEWTON:" << endl;

    newton_wrapper(target, max_number, x, regular_components, true);



    cout << "\nBIFURCATION ORDER:" << endl;

    int bound;

    vector<IVector> extra_special, new_special;

    for (intvec_it = special.begin(); intvec_it != special.end(); intvec_it++)
    {

        p = *intvec_it;

        if ( (bound = bifurcation_order(target, x, p, max_derivative, tolerance)) == -1 ) {

            cout << "Failed to find a bound for parameters " << p << ". Try a smaller tolerance or smaller parameter width." << endl;

            extra_special.push_back(p);

        } else if ( bound >= 0 ) {

            if (bound > max_number) extra_special.push_back(p);

            else new_special.push_back(p);

        }

    }

    special = new_special;

    cout << "Number of extra special boxes is " << extra_special.size() << endl;



    // write boxes to file streams

    int stream_precision = 6;

    if ( p.dimension() == 2 )
    {

        int reg_size = regular_components.size();

        ofstream regular_streams[reg_size];

        stringstream sstm;

        for (int i = 0; i < reg_size; i++)
        {

            sstm.str("");

            sstm << "components/regular" << i + 1;

            regular_streams[i].open(sstm.str());

            if (!regular_streams[i])
            {

                cerr << "cannot open file " << sstm.str() << '\n';
                exit(1);

            }

            regular_streams[i].setf(ios_base::fixed, ios_base::floatfield);

            regular_streams[i].precision(stream_precision);

            for (intvec_it = regular_components[i].begin(); intvec_it != regular_components[i].end(); intvec_it++)

                write_vector(*intvec_it, &(regular_streams[i]));

        }

        sstm.str("");
        sstm << "components/special";
        ofstream special_stream(sstm.str());
        if (!special_stream)
        {

            cerr << "cannot open file " << sstm.str() << '\n';
            exit(1);

        }
        special_stream.setf(ios_base::fixed, ios_base::floatfield);
        special_stream.precision(stream_precision);

        for (intvec_it = special.begin(); intvec_it != special.end(); intvec_it++) write_vector(*intvec_it, &special_stream);

        sstm.str("");
        sstm << "components/extra_special";
        ofstream extra_special_stream(sstm.str());
        if (!extra_special_stream)
        {

            cerr << "cannot open file " << sstm.str() << '\n';
            exit(1);

        }
        extra_special_stream.setf(ios_base::fixed, ios_base::floatfield);
        extra_special_stream.precision(stream_precision);

        for (intvec_it = extra_special.begin(); intvec_it != extra_special.end(); intvec_it++) write_vector(*intvec_it, &extra_special_stream);

        for (int i = 0; i < reg_size; i++) regular_streams[i].close();

        special_stream.close();

        extra_special_stream.close();

    }

}
*/


int automatic(IMap& target, IVector x, IVector p, int max_number, int max_derivative, double tolerance)
{

    State state;

    state.special.push_back(p);

    state.tolerance = tolerance;

    do {

        bisection(target, x, state); // subdivides special boxes

        bifurcation_order_wrapper(target, state, x, max_number, max_derivative);


        state.tolerance /= 2.;

    } while (!state.special.empty() && state.tolerance > 1e-5);

    if (!state.special.empty())

        return(AUTO_TOLERANCE);

    vector<vector<IVector>> regular_components;

    find_connected_components(state.regular, regular_components);

    return(newton_wrapper(target, regular_components, max_number, x, false));

}

int main()
{

    int max_derivative = 3; // maximum multiplicity

    int max_number = 5;

    IMap target = get_target(max_derivative);
    //IMap target("par:a,b,c;var:x;fun:(a*x)+(b*(x^2))+(c*(x^3));", max_derivative);

    IVector x(1), p(3);

    const double R = 0.3;
    interval radius = interval(-R, R);

    x[0] = interval(-R - 0.001, R);   // phase variable

    p[0] = interval(-1.) + radius;
    p[1] = interval(2.) + radius;
    p[2] = interval(1.) + radius;

    double tolerance = 1e-1;

    int result = automatic(target, x, p, max_number, max_derivative, tolerance); cout << "Result is " << result << endl;
    //manual(target, max_number, max_derivative, x, p, tolerance);

    return(0);

}
