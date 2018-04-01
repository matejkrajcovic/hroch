#include"random.h"
#include<random>
#include<unistd.h>

namespace likelihood {
default_random_engine generator;
uniform_real_distribution<double> distribution;

void random_init(){
    srand(time(NULL)+getpid());
    generator = default_random_engine(rand());
    distribution = uniform_real_distribution<double>();
}

double random_double(){
    return distribution(generator);
}
double random_double(double from, double to){
    return double(generator()-generator.min())/double(generator.max()-generator.min())*(to-from)+from;
}
}
