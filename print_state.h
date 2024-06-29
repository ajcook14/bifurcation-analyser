
#ifndef MANUAL_H
#define MANUAL_H

#include "capd/capdlib.h"
#include "solvers.h"
#include <iostream>
#include <fstream>

using namespace capd;
using namespace std;

void write_vector(IVector&, ofstream*);
void print_state(State& state);

#endif //MANUAL_H
