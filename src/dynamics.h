// Implementation of candidate sampling algorithm.

#ifndef DYNAMICS_H
#define DYNAMICS_H

#include"constants.h"
 

class Dynamics {
    const double startp = 1.0;
    History* history;
    HEvent* event;
    vector<int> atoms;
    vector<vector<vector<double>>> cherryness;
    vector<vector<vector<double>>> mass;
    vector<pair<double, trint>> all_endpoints;
    //vector<double> csum_ep;
    double sum_endpoints;
    double dupm, delm, deli, deli3;
    int n;
public:
    Dynamics(History* history, HEvent* event);
    void compute_graph(double dupm, double delm, double deli);
    Candidate get_candidate(bool maximal = false);
};

#endif
