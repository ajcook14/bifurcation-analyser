#include "capd/capdlib.h"

#include "./solvers.h"
#include "./neuron.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace capd;



void error(const char* p, const char* p2=" ")
{

    cerr << p << ' ' << p2 << '\n';
    exit(1);

}

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

    components.push_back(component);

}

int main()
{

    int max_derivative = 3; // maximum multiplicity

    int result;

    double tolerance;

    IMap target = get_target(max_derivative);

    //IMap target("par:a,b;var:x;fun:x^3-a*x+b;", max_derivative);

    IVector x(1), p(2);



    // initialise file streams

    ofstream regular_stream1("regular1");
    ofstream regular_stream2("regular2");
    ofstream special_stream("special");

    if (!regular_stream1) error("cannot open file", "regular1");
    if (!regular_stream2) error("cannot open file", "regular2");
    if (!special_stream) error("cannot open file", "special");

    regular_stream1.setf(ios_base::fixed, ios_base::floatfield);
    regular_stream2.setf(ios_base::fixed, ios_base::floatfield);
    special_stream.setf(ios_base::fixed, ios_base::floatfield);

    regular_stream1.precision(3);
    regular_stream2.precision(3);
    special_stream.precision(3);




    cout << "\nBISECTION:" << endl;

    vector<IVector> regular, special;

    tolerance = 1e-2;

    x[0] = interval(-1, 1);   // phase variable
    p[0] = interval(0, 2);
    p[1] = interval(-1, 1);

    cout << "\\x = " << x << " and \\p =  " << p << endl;
    bisection(target, x, p, regular, special, tolerance);

    cout << "Number of eliminated boxes = " << regular.size() << endl;
    cout << "Number of special boxes = " << special.size() << endl;

    vector<vector<IVector>> regular_components;

    find_connected_components(regular, regular_components);

    cout << "Number of components is " << regular_components.size() << endl;



    cout << "\nNEWTON:" << endl;

    int newton_verified[regular_components.size()]; // number of zeroes in each component

    x[0] = interval(-1, 1);

    for (int i = 0; i < regular_components.size(); i++)
    {

        p = regular_components[i][0];
        p = midVector(p);

        cout << "\\x = " << x << " and \\p = " << p << endl;

        while ( (result = newton_method(target, x, p)) == -1 )
        {

            cout << "Double counted solution" << endl;

            x[0] += interval(0, 1e-5);

            cout << "\\x = " << x << " and \\p = " << p << endl;

        }

        if ( result == -2 ) cout << "Couldn't resolve solutions - try a smaller inflation or different parameters" << endl;

        else {

            cout << "Newton-verified number of solutions in component " << i << " is " << result << "." << endl;

            newton_verified[i] = result;

        }

    }

    vector<IVector>::iterator intvec_it;

    for (intvec_it = regular_components[0].begin(); intvec_it != regular_components[0].end(); intvec_it++) write_vector(*intvec_it, &regular_stream1);
    for (intvec_it = regular_components[1].begin(); intvec_it != regular_components[1].end(); intvec_it++) write_vector(*intvec_it, &regular_stream2);
    for (intvec_it = special.begin(); intvec_it != special.end(); intvec_it++) write_vector(*intvec_it, &special_stream);




    cout << "\nBIFURCATION ORDER:" << endl;

    x[0] = interval(-1, 1);
    p[0] = interval(1-1e-2, 1+1e-2);
    p[1] = interval(1e-5, 2e-5);

    cout << "\\x =  " << x << " and \\p = " << p << endl;

    tolerance = 1e-1;

    int bound;

    if ( (bound = bifurcation_order(target, x, p, max_derivative, tolerance)) == -1 ) {

        cout << "Failed to find a bound. Try a smaller tolerance or smaller parameter width." << endl;

    } else if ( bound >= 0 ) {

        cout << "Upper bound on number of zeroes is " << bound << endl;

    }



    regular_stream1.close();
    regular_stream2.close();
    special_stream.close();

    return(0);

}
