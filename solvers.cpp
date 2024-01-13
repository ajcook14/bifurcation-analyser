#include "capd/capdlib.h"
#include "capd/newton/Newton.h"

#include "small_funcs.h"

#include <iostream>
#include <queue>

using namespace std;
using namespace capd;

void cartesian(IVector& a, IVector& b, IVector& product)
{

    /**
        Concatenates a and b and puts the result in product.
    **/

    IVector temp(a.dimension() + b.dimension());

    IVector::iterator intvec_it;

    for (int i = 0; i < a.dimension(); i++) temp[i] = a[i];

    for (int i = 0; i < b.dimension(); i++) temp[i + a.dimension()] = b[i];

    product = temp;

}



int bisection_aux(IMap& target, IVector x, IVector& p, double tolerance)
{

    /**
        Returns 0 if the region x in phase space does not contain any bifurcations for all parameters in p.
        Returns 1 if the region p could not be eliminated (possibly because there is a bifurcation).
    **/

    IVector left, right;

    for (int i = 0; i < p.dimension(); i++) target.setParameter(i, p[i]);

    queue<IVector> bisection_queue;

    bisection_queue.push(x);

    // interval bisection

    while (!bisection_queue.empty())
    {

        x = bisection_queue.front();

        bisection_queue.pop();

        if ( !containsZero(target(x)) || !containsZero(target.derivative(x)) ) continue;

        if (maxDiam(x) < tolerance) return(1);

        else {

            subdivide(x, left, right);

            bisection_queue.push(left);
            bisection_queue.push(right);

        }

    }

    return(0);

}


void bisection(IMap& target, IVector x, IVector p, vector<IVector>& regular, vector<IVector>& special, double tolerance)
{

    /**
        Eliminates all intervals in parameter space larger than tolerance that don't contain bifurcations and places them in regular.
        Places remaining intervals in special.
    **/

    IVector left, right;

    queue<IVector> bisection_queue;

    bisection_queue.push(p);

    // interval bisection

    while (!bisection_queue.empty())
    {

        p = bisection_queue.front();

        bisection_queue.pop();

        if ( bisection_aux(target, x, p, tolerance) == 0 ) {

            regular.push_back(p);

            continue;

        }

        if (maxDiam(p) < tolerance) special.push_back(p);

        else {

            subdivide(p, left, right);

            bisection_queue.push(left);
            bisection_queue.push(right);

        }

    }

}



int newton_method(IMap& target, IVector& x, IVector& p)
{

    /**
        Returns the exact number of solutions in the initial interval x or a negative value if solutions could not be resolved.
    **/

    for (int i = 0; i < p.dimension(); i++) target.setParameter(i, p[i]);

    IVector midpoint, N, left, right, epsilon(1);

    epsilon[0] = interval(-1e-10, 1e-10);

    queue<IVector> bisection_queue;

    bisection_queue.push(x);

    // interval bisection

    list<IVector> verified;

    while ( !bisection_queue.empty() )
    {

        x = bisection_queue.front();

        bisection_queue.pop();

        if ( min_dim(x) < 4 * width(epsilon[0]) ) return(-2);   // check min dimension of x is at least twice width(epsilon)

        if ( !containsZero(target(x)) ) continue;

        if ( containsZero(target.derivative(x)) ) {

            subdivide(x, left, right);

            bisection_queue.push(left);
            bisection_queue.push(right);

            continue;

        }

        midpoint = midVector(x);
        N = newton::NewtonOperator(midpoint, x, target);

        if ( subset(N, x) ) {

            verified.push_back(x);
            continue;

        }

        if ( intersectionIsEmpty(N, x) ) continue;

        // all other cases:

        subdivide(x, left, right);

        left += epsilon;
        right += epsilon;

        bisection_queue.push(left);
        bisection_queue.push(right);

    }

    // check verified list is pairwise disjoint
    list<IVector>::iterator i, j;

    for (i = verified.begin(); i != verified.end(); i++) {

        j = i;
        j++;

        for (; j != verified.end(); j++) {

            if ( intersectionIsEmpty(*i, *j) ) continue;

            if ( !containsZero(target(intersection(*i, *j))) ) continue;

            return(-1); // double counted solution

        }

    }

    return(verified.size());

}



double total_measure(vector<IVector>& intervals)    // assumes the measure of the union equals the sum of the measures
{

    double measure = 0;

    vector<IVector>::iterator i;

    for (i = intervals.begin(); i != intervals.end(); i++) {

        measure += diam(*i)[0].leftBound();

    }

    return(measure);

}



void refine_measure(IMap& target, vector<IVector>& intervals, double tolerance)
{

    /**
        Subdivides all intervals and throws away any infeasible intervals.
        Halts when the change in total measure is smaller than tolerance.
        The tolerance should be chosen so it only halts when the total measure
        plateaus due to the fact that target is a parameterised family of functions.
    **/

    double measure, difference;

    vector<IVector> temp;

    vector<IVector>::iterator i;

    IVector left(1), right(1);

    do {

        measure = total_measure(intervals);

        temp.clear();

        for (i = intervals.begin(); i != intervals.end(); i++) {

            subdivide(*i, left, right);

            temp.push_back(left);
            temp.push_back(right);

        }

        intervals.clear();

        for (i = temp.begin(); i != temp.end(); i++) {

            if ( containsZero(target(*i)) ) intervals.push_back(*i);

        }

        difference = measure - total_measure(intervals);

    } while ( difference > tolerance || difference == 0 );

}



void connected_components(vector<IVector>& intervals, vector<vector<IVector>>& components)
{

    /**
        Assumes the IVectors in intervals are dimension 1 only.
        Writes the connected components of intervals to components.
        The vector intervals will get sorted into connected components.
    **/

    std::sort(intervals.begin(), intervals.end(), [](IVector &a, IVector &b) { return(midVector(a)[0].leftBound() < midVector(b)[0].leftBound()); });

    vector<IVector> component;

    component.push_back(intervals[0]);

    vector<IVector>::iterator interval_it;

    for (int i = 0; i < intervals.size() - 1; i++)
    {

        if ( !intersectionIsEmpty( intervals[i], intervals[i + 1] ) ) component.push_back(intervals[i + 1]);//component = intervalHull(component, intervals[i + 1]);

        else {

            components.push_back(component);

            component.clear();

            component.push_back(intervals[i + 1]);

        }

    }

    components.push_back(component);    // append remaining component

}



int bifurcation_order(IMap& target, IVector& x, IVector& p, int max_derivative, double tolerance)
{

    /**
        Returns an upper bound for the number of solutions to target(x) = 0 over the region p of parameters.
        Requires knowledge of max_derivative, the maximum order of each solution.
        tolerance should be of similar order of magnitude to the diameter of p. See documentation for refine_measure.
    **/

    for (int i = 0; i < p.dimension(); i++) target.setParameter(i, p[i]);

    vector<IVector> feasible;

    feasible.push_back(x);

    refine_measure(target, feasible, tolerance);

    vector<vector<IVector>> components;

    connected_components(feasible, components);



    int dimOut, dimIn;
    dimOut = dimIn = 1;

    IJet jet(dimOut, dimIn, max_derivative);

    Multipointer mp;

    vector<vector<IVector>>::iterator comp_it;
    vector<IVector>::iterator interval_it;

    int estimate = 0;

    int n_intervals = 0;

    for (comp_it = components.begin(); comp_it != components.end(); comp_it++)
    {

        for (int order = 1; order <= max_derivative + 1; order++)
        {

            if ( order == max_derivative + 1 ) return(-1);

            n_intervals = 0;

            for (interval_it = (*comp_it).begin(); interval_it != (*comp_it).end(); interval_it++)
            {

                target(*interval_it, jet);

                mp = jet.first(order);

                if ( !containsZero(jet(mp)) ) n_intervals++;

            }

            if ( n_intervals == (*comp_it).size() ) {

                estimate += order;

                break;

            }

        }

    }

    return(estimate);

}

/*
int bifurcation_order(IMap& target, IVector& x, IVector& p, vector<IVector>& feasible, int maxDerivative, double TOL)
// returns an upper bound on the number of solutions of target in box
{

    IVector left(1), right(1);

    int dimOut, dimIn;
    dimOut = dimIn = 1;

    IJet jet(dimOut, dimIn, maxDerivative);

    Multipointer mp;

    target.setParameter(0, p[0]);


    // initialise queue

    queue<IVector> bisection_queue;

    bisection_queue.push(x);

    // the following loop produces a cpp vector of intervals (feasible) whose union contains all zeroes of the target, and such that
    // at least one derivative of order < maxDerivative is non-vanishing on each interval.

    int order;

    while (!bisection_queue.empty())
    {

        x = bisection_queue.front();

        bisection_queue.pop();

        if ( !containsZero(target(x)) ) continue;

        target(x, jet);

        for (order = 0; order <= maxDerivative; order++)
        {

            mp = jet.first(order);
            if ( !containsZero(jet(mp)) ) {

                feasible.push_back(x);

                break;

            }


            if ( order == maxDerivative ) {

                if ( maxDiam(x) < TOL ) return(-1);

                subdivide(x, left, right);

                bisection_queue.push(left);
                bisection_queue.push(right);

            }

        }

    }

    // check for connected components

    std::sort(feasible.begin(), feasible.end(), [](IVector &a, IVector &b) { return(midVector(a)[0].leftBound() < midVector(b)[0].leftBound()); });

    int components[feasible.size()];

    int ci = 0; //component index

    components[0] = 0;

    for (int i = 0; i < feasible.size() - 1; i++)
    {

        if ( intersectionIsEmpty( feasible[i], feasible[i + 1] ) ) ci++;

        components[i + 1] = ci;

    }

    int c_num = ci + 1; // number of componenets

    cout << "c_num = " << c_num << endl;

    std::for_each(feasible.begin(), feasible.end(), [](IVector const v) { cout << v << endl; });

    int orders[c_num];

    bool break_order;

    int ii = 0; // interval index

    int ii_begin = 0;

    cout << "feasible.size = " << feasible.size() << endl;

    for (int i = 0; i < feasible.size(); i++) cout << components[i];
    cout << "\n";

    for (ci = 0; ci < c_num; ci++)
    {

        break_order = false;

        ii_begin = ii;

        for (order = 1; order <= maxDerivative; order++)
        {

            ii = ii_begin;

            target(feasible[ii], jet);
            mp = jet.first(order);

            while ( !containsZero(jet(mp)) )
            {

                //cout << "ii = " << ii << endl;

                ii++;

                if ( ii >= feasible.size() || components[ii] != ci ) {

                    orders[ci] = order;

                    break_order = true;

                    break;

                }

                target(feasible[ii], jet);
                mp = jet.first(order);

            }

            if ( break_order ) break;

        }

        if ( !break_order ) {

            cout << "bifurcation_order: max derivative exceeded" << endl;

            return(-1);

        }

    }


    for (int i = 0; i < c_num; i++) {

        cout << "component = " << i << ", order = " << orders[i] << endl;

    }

    return(0);

}
*/
