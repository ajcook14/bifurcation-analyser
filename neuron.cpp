#include "capd/capdlib.h"
#include "neuron.h"

using namespace capd;

using capd::autodiff::Node;



Node tanh(Node x) {

    return 1 - (2/(exp(2*x) + 1));

}

Node cosh(Node x) {

    return (exp(x) + exp(-x))/2;

}

Node dtanh(Node x) {

    return( 1/((cosh(x))^2) );

}

void neuron(Node t,
            Node in[], int dimIn,
            Node out[], int dimOut,
            Node params[], int noParam)
{

    assert(dimIn == 1);
    assert(dimOut == 1);
    assert(noParam == 2);

    Node x = in[0];
    Node a = params[0];
    Node b = params[1];

    //out[0] = tanh((a * x) + b) - x;
    //out[0] = (a * dtanh(a * (x - b))) - (1./2.);
    out[0] = tanh(2 * (x - a)) + tanh(b * (x - (-5.))) - x;
    //out[0] = (a * x) + (b * (x^2)) + (-0.9375 * (x^3));

}

void small_net(Node t,
               Node in[], int dimIn,
               Node out[], int dimOut,
               Node params[], int noParam)
{

    assert(dimIn == 1);
    assert(dimOut == 1);
    assert(noParam == 3);

    Node& x = in[0];
    Node& a = params[0];
    Node& b = params[1];
    Node& c = params[2];

    //out[0] = (a * dtanh((a * x) + b)) + (c * dtanh(c * x)) - 1;
    //out[0] = (a * x) + (b * (x^2)) + (c * (x^3));
    out[0] = tanh(2 * (x - a)) + tanh(b * (x - c)) - x;

}

IMap get_target(int maxDerivative)
{

    int dimIn, dimOut, noParam;
    dimIn = dimOut = 1;
    noParam = 3;

    IMap target(small_net, dimIn, dimOut, noParam, maxDerivative);

    return(target);

}