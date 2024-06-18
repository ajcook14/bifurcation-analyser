
#ifndef MANUAL_H
#define MANUAL_H

using namespace capd;
using namespace std;

void write_vector(IVector&, ofstream*);
void manual(IMap& target, int max_number, int max_derivative, IVector x, IVector p, double tolerance);

#endif //MANUAL_H
