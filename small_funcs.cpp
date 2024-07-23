#include "capd/capdlib.h"
#include "small_funcs.h"

using namespace std;
using namespace capd;



double min_dim(IVector& x)
{
    double min_dim = widths(x)[0];
    for (std::size_t i = 0; i < x.dimension(); i++) {
        if ( widths(x)[i] < min_dim ) {
            min_dim = widths(x)[i];
        }
    }

    return(min_dim);
}

unsigned int arg_max_dim(IVector& x)
{
    unsigned int index_max_dim = 0;

    for (std::size_t i = 0; i < x.dimension(); i++) {
        if ( widths(x)[i] > widths(x)[index_max_dim] ) {
            index_max_dim = i;
        }
    }
    return(index_max_dim);
}

void subdivide(IVector& x, IVector& left, IVector& right)
{
    left = x;
    right = x;

    unsigned int max_dim = arg_max_dim(x);

    left[max_dim] = interval(x[max_dim].leftBound(), mid(x[max_dim]).rightBound());
    right[max_dim] = interval(mid(x[max_dim]).leftBound(), x[max_dim].rightBound());
}
