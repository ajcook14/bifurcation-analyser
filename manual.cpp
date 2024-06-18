
#include "capd/capdlib.h"


#include <iostream>
#include <fstream>
#include "manual.h"
#include "solvers.h"

using namespace std;

void write_vector(IVector& iv, ofstream* to)
{
    /**
        assumes iv is dimension 2
    **/

    *to << iv[0].leftBound() << ',' << iv[0].rightBound() << ','
       << iv[1].leftBound() << ',' << iv[1].rightBound() << endl;
}

void manual(IMap& target, int max_number, int max_derivative, IVector x, IVector p, double tolerance)
{
    State state;
    state.special.push_back(p);
    state.tolerance = tolerance;

    cout << "\nBISECTION:" << endl;
    vector<IVector> regular, special;
    cout << "\\x = " << x << " and \\p =  " << p << endl;
    bisection(target, x, p, regular, special, tolerance);
    cout << "Number of eliminated boxes = " << regular.size() << endl;
    cout << "Number of special boxes = " << special.size() << endl;
    vector<vector<IVector>> regular_components;
    find_connected_components(state.regular, regular_components);
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