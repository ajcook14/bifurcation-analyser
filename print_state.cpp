
#include "capd/capdlib.h"


#include <iostream>
#include <fstream>
#include "print_state.h"
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

void print_state(State& state)
{
    // write boxes to file streams
    int stream_precision = 6;
    if (!state.special.empty()) {
        assert(state.special[0].dimension() == 2);
    }

    vector<vector<IVector>> regular_components;
    find_connected_components(state.regular, regular_components);
    std::size_t reg_size = regular_components.size();
    ofstream regular_streams[reg_size];
    stringstream sstm;
    for (std::size_t i = 0; i < reg_size; i++)
    {
        sstm.str("");
        sstm << "../components/regular" << i + 1;
        regular_streams[i].open(sstm.str());
        if (!regular_streams[i]) {
            cerr << "cannot open file " << sstm.str() << '\n';
            exit(1);
        }
        regular_streams[i].setf(ios_base::fixed, ios_base::floatfield);
        regular_streams[i].precision(stream_precision);
        for (auto & iv : regular_components[i]) {
            write_vector(iv, &(regular_streams[i]));
        }
    }
    sstm.str("");
    sstm << "../components/verified";
    ofstream verified_stream(sstm.str());
    if (!verified_stream)
    {
        cerr << "cannot open file " << sstm.str() << '\n';
        exit(1);
    }
    verified_stream.setf(ios_base::fixed, ios_base::floatfield);
    verified_stream.precision(stream_precision);
    for (auto & iv : state.verified) {
        write_vector(iv, &verified_stream);
    }
    sstm.str("");
    sstm << "../components/special";
    ofstream special_stream(sstm.str());
    if (!special_stream)
    {
        cerr << "cannot open file " << sstm.str() << '\n';
        exit(1);
    }
    special_stream.setf(ios_base::fixed, ios_base::floatfield);
    special_stream.precision(stream_precision);
    for (auto & iv : state.special) {
        write_vector(iv, &special_stream);
    }
    for (std::size_t i = 0; i < reg_size; i++) {
        regular_streams[i].close();
    }
    verified_stream.close();
    special_stream.close();
}