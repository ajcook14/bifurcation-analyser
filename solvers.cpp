#include "capd/capdlib.h"
#include "capd/newton/Newton.h"

#include "small_funcs.h"
#include "solvers.h"
#include "capd_error.h"

#include <iostream>
#include <queue>

using namespace std;
using namespace capd;



int bisection_aux(IMap& target, IVector x, IVector p, double tolerance)
{

    /**
        Returns 0 if the region x in phase space does not contain any bifurcations for all parameters in p.
        Returns 1 if the region p could not be eliminated (possibly because there is a bifurcation).
    **/

    IVector left, right;

    for (int i = 0; i < p.dimension(); i++)

        target.setParameter(i, p[i]);

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


void bisection(IMap& target, IVector x, State& state)
{

    /**
        Subdivides intervals in state.special and eliminates all intervals larger than state.tolerance that don't contain bifurcations and places them in state.regular.
        Places remaining intervals in state.special.
    **/

    IVector left, right, p;

    vector<IVector> new_special;

    queue<IVector> bisection_queue;

    for_each (state.special.begin(), state.special.end(), [&bisection_queue](IVector p) {

        bisection_queue.push(p);

    });

    // interval bisection

    while (!bisection_queue.empty())
    {

        p = bisection_queue.front();

        bisection_queue.pop();

        if ( bisection_aux(target, x, p, state.tolerance) == 0 ) {

            state.regular.push_back(p);

        } else if (maxDiam(p) < state.tolerance) {

            new_special.push_back(p);

        } else {

            subdivide(p, left, right);

            bisection_queue.push(left);
            bisection_queue.push(right);

        }

    }

    state.special.clear();

    state.special = new_special;

}



int newton_method(IMap& target, IVector x, IVector p)
{

    /**
        Returns the exact number of solutions in the initial interval x or a negative value if solutions could not be resolved.
        Assumes there are no bifurcations in p.
    **/

    for (int i = 0; i < p.dimension(); i++) target.setParameter(i, p[i]);

    IVector midpoint, N, left, right, epsilon(1);

    epsilon[0] = interval(-1e-17, 1e-17);

    queue<IVector> bisection_queue;

    bisection_queue.push(x);

    // interval bisection

    list<IVector> verified;

    while ( !bisection_queue.empty() )
    {

        x = bisection_queue.front();

        bisection_queue.pop();

        if ( min_dim(x) < 4 * width(epsilon[0]) ) return(NEWTON_EPSILON_WIDTH);   // check min dimension of x is at least twice width(epsilon)

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

        if ( subset(x, N) ) {

            subdivide(x, left, right);

            bisection_queue.push(left);
            bisection_queue.push(right);

        }

        //left += epsilon;
        //right += epsilon;

        bisection_queue.push(intersection(N, x));

    }

    // check verified list is pairwise disjoint
    list<IVector>::iterator i, j;

    for (i = verified.begin(); i != verified.end(); i++) {

        j = i;
        j++;

        for (; j != verified.end(); j++) {

            if ( intersectionIsEmpty(*i, *j) ) continue;

            if ( !containsZero(target(intersection(*i, *j))) ) continue;

            return(NEWTON_DOUBLE_COUNT); // double counted solution

        }

    }

    return(verified.size());

}



double total_measure(vector<IVector>& intervals)    // assumes intervals are pairwise disjoint
{

    double measure = 0;

    vector<IVector>::iterator i;

    for (i = intervals.begin(); i != intervals.end(); i++) {

        measure += diam(*i)[0].leftBound();

    }

    return(measure);

}



int refine_measure(IMap& target, vector<IVector>& intervals, double tolerance)
{

    /**
        Subdivides all intervals and throws away any infeasible intervals.
        Halts when the change in total measure is smaller than tolerance, or
        when there are no intervals left, or when no intervals have been removed
        after REFINE_MAX_SUBDIVISIONS subdivisions.
        The tolerance should be chosen so it only halts when the total measure
        plateaus.
    **/

    double measure, current_measure, difference;

    vector<IVector> temp;

    vector<IVector>::iterator i;

    IVector left(1), right(1);

    double init_measure = total_measure(intervals);

    int timeout = 0;

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

        if (intervals.empty())

            return(ERROR_NO_BIFURCATION);

        current_measure = total_measure(intervals);

        if (timeout == REFINE_MAX_SUBDIVISIONS && init_measure - current_measure == 0.)

            return(ERROR_MAX_SUBDIVISIONS);

        timeout++;

        difference = measure - current_measure;

    } while ( difference > tolerance || difference == 0 );

    return(0);

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





int bifurcation_order(IMap& target, IVector x, IVector p, int max_derivative, double tolerance)
{

    /**
        Returns an upper bound for the number of solutions to target(x) = 0 over the region p of parameters.
        Requires knowledge of max_derivative, the maximum order of each solution.
        tolerance should be of similar order of magnitude to the diameter of p. See documentation for refine_measure.
    **/

    for (int i = 0; i < p.dimension(); i++) target.setParameter(i, p[i]);

    vector<IVector> feasible;

    feasible.push_back(x);

    int error;

    if ( (error = refine_measure(target, feasible, tolerance)) < 0 )

        return(error);

    vector<vector<IVector>> components;

    find_connected_components(feasible, components);



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

            if ( order == max_derivative + 1 )

                return(ERROR_MAX_DERIVATIVE);

            n_intervals = 0;

            for (interval_it = comp_it->begin(); interval_it != comp_it->end(); interval_it++)
            {

                target(*interval_it, jet);

                mp = jet.first(order);

                if ( !containsZero(jet(mp)) ) n_intervals++;

            }

            if ( n_intervals == comp_it->size() ) {

                estimate += order;

                break;

            }

        }

    }

    return(estimate);

}

