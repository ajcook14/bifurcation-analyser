
#ifndef BIFURCATION_ANALYSER_SMALL_FUNCS_H
#define BIFURCATION_ANALYSER_SMALL_FUNCS_H

#endif //BIFURCATION_ANALYSER_SMALL_FUNCS_H

using namespace std;
using namespace capd;



double min_dim(IVector&);

int arg_max_dim(IVector&);

void subdivide(IVector& x, IVector& left, IVector& right);

void subdivide(SpecialBox& s_box, SpecialBox& left, SpecialBox& right);