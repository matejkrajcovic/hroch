// Implementation of candidate sampling algorithm.

#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <vector>
#include <utility>
#include <tuple>

#include "history.h"
#include "hevent.h"

class Dynamics {
    const double startp = 1.0;
    History* history;
    HEvent* event;
    std::vector<int> atoms;
    std::vector<std::vector<std::vector<double>>> cherryness;
    std::vector<std::vector<std::vector<double>>> mass;
    std::vector<std::pair<double, std::tuple<int, int, int>>> all_endpoints;
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
