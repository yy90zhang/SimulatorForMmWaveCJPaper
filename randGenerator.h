#ifndef RANDGENERATOR_H_INCLUDED
#define RANDGENERATOR_H_INCLUDED

#include <random>
using namespace std;



///generating uniform rv between a and b
double uniform_rv_generator(double a, double b);

///generating poisson rv with mean mean
int poisson_rv_generator(double mean);

///generating gamma rv with shape alpha and scale beta (i.e., rate 1/beta)
double gamma_rv_generator(double alpha, double beta);

#endif // RANDGENERATOR_H_INCLUDED
