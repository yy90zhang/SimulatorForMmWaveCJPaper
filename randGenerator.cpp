#include "randGenerator.h"

random_device rd;
default_random_engine seed(rd());
/*Note these functions are Pseudo-randomÅC which means that the
simulation results (e.g., probabilities) of events will not change with time*/

///generating uniform rv between a and b
double uniform_rv_generator(double a, double b)
{
    uniform_real_distribution<double> distribution(a, b);
    return distribution(seed);
}

int poisson_rv_generator(double mean)
{
    poisson_distribution<int> distribution(mean);
    return distribution(seed);
}

///generating gamma rv with shape alpha and scale beta (i.e., rate 1/beta)
double gamma_rv_generator(double alpha, double beta)
{
    gamma_distribution<double> distribution(alpha, beta);
    return distribution(seed);
}
