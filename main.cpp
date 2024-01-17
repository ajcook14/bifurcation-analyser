#include "capd/capdlib.h"

#include "./solvers.h"
#include "./neuron.h"

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

    components.push_back(component);

}

int main()
{

    int max_derivative = 3; // maximum multiplicity

    int result;

    double tolerance;

    IMap target = get_target(max_derivative);

    //IMap target("par:a,b;var:x;fun:x^3-a*x+b;", max_derivative);
    //IMap target("par:a,b;var:x;fun:x^4+a*(x^2)+b*x;", max_derivative);

    IVector x(1), p(2);

    vector<IVector>::iterator intvec_it;







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

    IVector x_newton;

    for (int i = 0; i < regular_components.size(); i++)
    {

        x_newton = x;

        p = regular_components[i][0];
        p = midVector(p);

        cout << "\\x = " << x_newton << " and \\p = " << p << endl;

        while ( (result = newton_method(target, x_newton, p)) == -1 )
        {

            cout << "Double counted solution" << endl;

            x_newton[0] += interval(0, 1e-5);

            cout << "\\x = " << x_newton << " and \\p = " << p << endl;

        }

        if ( result == -2 ) cout << "Couldn't resolve solutions - try a smaller inflation or different parameters" << endl;

        else {

            cout << "Newton-verified number of solutions in component " << i << " is " << result << "." << endl;

            newton_verified[i] = result;

        }

    }



    cout << "\nBIFURCATION ORDER:" << endl;

    tolerance = 1e-2;

    int bound;

    vector<IVector> extra_special, new_special;

    for (intvec_it = special.begin(); intvec_it != special.end(); intvec_it++)
    {

        p = *intvec_it;

        if ( (bound = bifurcation_order(target, x, p, max_derivative, tolerance)) == -1 ) {

            cout << "Failed to find a bound for parameters " << p << ". Try a smaller tolerance or smaller parameter width." << endl;

            extra_special.push_back(p);

        } else if ( bound >= 0 ) {

            if (bound > max_derivative) extra_special.push_back(p);

            else new_special.push_back(p);

        }

    }

    special = new_special;

    cout << "Number of extra special boxes is " << extra_special.size() << endl;

    //cout << "Upper bound on number of zeroes is " << max_bound << endl;



    // write boxes to file streams
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

            regular_streams[i].precision(3);

            for (intvec_it = regular_components[i].begin(); intvec_it != regular_components[i].end(); intvec_it++) write_vector(*intvec_it, &(regular_streams[i]));

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
        special_stream.precision(3);

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
        extra_special_stream.precision(3);

        for (intvec_it = extra_special.begin(); intvec_it != extra_special.end(); intvec_it++) write_vector(*intvec_it, &extra_special_stream);

        for (int i = 0; i < reg_size; i++) regular_streams[i].close();

        special_stream.close();

        extra_special_stream.close();

    }

    return(0);

}
