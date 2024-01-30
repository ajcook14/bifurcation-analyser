#include "capd/capdlib.h"

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

    Node x = in[0];
    Node a = params[0];
    Node b = params[1];

    //out[0] = tanh((a * x) + b) - x;
    //out[0] = (a * dtanh(a * (x - b))) - (1./2.);
    out[0] = tanh(2 * (x + 1)) + tanh(a * (x - b)) - x;

}

IMap get_target(int maxDerivative)
{

    int dimIn, dimOut, noParam;
    dimIn = dimOut = 1;
    noParam = 2;

    //IMap target("par:a;var:x;fun:a*x-x^3;", maxDerivative);

    IMap target(neuron, dimIn, dimOut, noParam, maxDerivative);

    return(target);

}
