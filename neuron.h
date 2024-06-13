
#ifndef BIFURCATION_ANALYSER_NEURON_H
#define BIFURCATION_ANALYSER_NEURON_H

#endif //BIFURCATION_ANALYSER_NEURON_H

using capd::autodiff::Node;
using namespace capd;

void neuron(Node t, Node in[], int dimIn, Node out[], int dimOut, Node params[], int noParam);

IMap get_target(int maxDerivative);