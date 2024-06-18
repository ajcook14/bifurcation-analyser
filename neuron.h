#ifndef NEURON_H
#define NEURON_H
using capd::autodiff::Node;
using namespace capd;

void neuron(Node t, Node in[], int dimIn, Node out[], int dimOut, Node params[], int noParam);

IMap get_target(int maxDerivative, int phase_dim, int param_dim);
#endif // NEURON_H
