#include "capd/capdlib.h"
#include "capd/newton/Newton.h"

#include "solvers.h"
#include "small_funcs.h"
#include "capd_error.h"

#include <iostream>
#include <queue>
#include <iomanip>

using namespace std;
using namespace capd;



extern Statistics statistics;

void bisection(IMap&/******/ target, State& state)
{

    /**
        Subdivides intervals in state.special and eliminates all intervals larger than state.tolerance that don't contain bifurcations and places them in state.regular.
        Places remaining intervals in state.special.
    **/

    SpecialBox s_box, left, right;
    vector<SpecialBox> new_special;
    queue<SpecialBox> bisection_queue;

    for_each (state.special.begin(), state.special.end(), [&](SpecialBox s_box) {
        bisection_queue.push(s_box);
    });

    int result;
    while (!bisection_queue.empty()) {
        std::cout << '\r';
        std::cout << std::setw(10) << std::right << state.regular.size()
                  << std::setw(10) << std::right << state.special.size()
                  << std::setw(10) << std::right << state.verified.size()
                  << std::setw(13) << std::right << new_special.size() << "  "
                  << std::setw(15) << std::left << state.tolerance << std::flush;

        s_box = bisection_queue.front();
        bisection_queue.pop();


        result = refine_measure(target, &s_box, state.tolerance);

        if (result == ERROR_NO_BIFURCATION) {
            statistics.no_bifurcation++;
        } else if ( result == REACHED_MAX_SUBDIVISIONS ) {
            statistics.max_subdivisions++;
        }

        if ( result == ERROR_NO_BIFURCATION ) {
            state.regular.push_back(s_box.parameters);
        } else if ( maxDiam(s_box.parameters) < state.tolerance ) {
            new_special.push_back(s_box);
        } else {
            subdivide(s_box, left, right);
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

    // check verified list is pairwise disjoint/////////////////////put this into a separate subroutine
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

int refine_measure(IMap&/******/ target, SpecialBox* s_box, double tolerance) ////rename to divide_phase
{

    /**
        Subdivides s_box->phase_boxes and throws away any infeasible intervals.
        Halts when the change in total measure is smaller than tolerance, or
        when there are no intervals left, or when no intervals have been removed
        after REFINE_MAX_SUBDIVISIONS subdivisions.
        The tolerance should be chosen so it only halts when the total measure
        plateaus.
    **/

    for (int i = 0; i < s_box->parameters.dimension(); i++) target.setParameter(i, s_box->parameters[i]);

    double measure, current_measure, difference;

    vector<IVector> temp;

    IVector left, right;

    double init_measure = total_measure(s_box->phase_boxes);

    int timeout = 0;

    do {

        measure = total_measure(s_box->phase_boxes);

        temp.clear();

        for_each (s_box->phase_boxes.begin(), s_box->phase_boxes.end(), [&](IVector iv) {

            subdivide(iv, left, right);

            temp.push_back(left);
            temp.push_back(right);

        });

        s_box->phase_boxes.clear();

        for_each (temp.begin(), temp.end(), [&](IVector iv) {
            if ( containsZero(target(iv)) )
                s_box->phase_boxes.push_back(iv);
        });

        if (s_box->phase_boxes.empty())
            return(ERROR_NO_BIFURCATION);/////rename to NO_BIFURCATION

        current_measure = total_measure(s_box->phase_boxes);

        if (timeout == REFINE_MAX_SUBDIVISIONS && init_measure - current_measure == 0.)
            return(REACHED_MAX_SUBDIVISIONS);

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

int bifurcation_order(IMap&/******/ target, SpecialBox s_box, int max_derivative)
{

    /**
        Returns an upper bound for the number of solutions to target(x) = 0 over the s_box.
        Requires knowledge of max_derivative, the maximum order of each solution.
    **/

    for (int i = 0; i < s_box.parameters.dimension(); i++) target.setParameter(i, s_box.parameters[i]);

    vector<vector<IVector>> components;
    find_connected_components(s_box.phase_boxes, components);

    int dimOut, dimIn;
    dimOut = dimIn = 1;
    IJet jet(dimOut, dimIn, max_derivative);
    Multipointer mp;

    int estimate = 0;
    int n_intervals = 0;
    bool max_flag = false;

    for_each (components.begin(), components.end(), [&](vector<IVector> comp) {

        for (int order = 1; order <= max_derivative + 1; order++)
        {
            if ( order == max_derivative + 1 ) {
                max_flag = true;
                return;
            }

            n_intervals = 0;

            for_each (comp.begin(), comp.end(), [&](IVector iv) {
                target(iv, jet);
                mp = jet.first(order);
                if ( !containsZero(jet(mp)) ) n_intervals++;
            });

            if ( n_intervals == comp.size() ) {
                estimate += order;
                break;
            }
        }
    });

    if (max_flag)
        return(ERROR_MAX_DERIVATIVE);
    return(estimate);
}
