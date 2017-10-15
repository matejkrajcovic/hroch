// Probabilistic model used for generated data

#ifndef RANDOM_H
#define RANDOM_H

#include"constants.h"

void random_init();
double random_double(); // random from interval [0.0, 1.0)
double random_double(double from, double to); // random from interval [from, to)
int random_int(int from, int to); // random integer from interval [from, to)

class Model {
private:
    geometric_distribution<int> length_distribution;
    geometric_distribution<int> distance_distribution;
    exponential_distribution<double> time_distribution;

    const double prob_del = 0.05;
    const double prob_inv = 0.3902;
    const int mean_len = 14307;
    const int mean_dist = 306718;
    const double event_rate = 200;
    const double mut_alpha = 1./3.;
   
    static Model* _instance;
    Model();
public:
    static Model* instance();
    
    const int length_threshold = 50;
    bool is_random_del();
    bool is_random_inv();
    int get_random_len();
    int get_random_dist();
    double get_random_time();
    char get_mutated_base(char base, double time);
    bool get_indel_happened(double time);

    GEvent* get_random_event(int sequence_length);
};

#endif
