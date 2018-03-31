#include <string>
#include <fstream>
#include <cmath>
#include <cassert>
#include <boost/functional/hash.hpp>
#include <utility>

#include "machine.h"
#include "score.h"

using namespace std;

double MachineBachelor::predict(const vector<double>& values) {
    double sum =
        values[SC_DEL_NUM] * -1. +
        values[SC_EV_LEN] * 1. +
        values[SC_EV_SIDES] * -2. +
        values[SC_EV_DSIG] * 6.0 +
        (values[SC_EV_POST_BP] - values[SC_EV_PREV_BP]) * 6.0;
    return exp(sum);
}


void MachineLinear::train_data(const vector<double>& values, double result) {
    data.push_back(values);
    data.back().push_back(result);
}
void MachineLinear::save() {
    string strict = "";
    if (strict_compare) strict = "-strict";
    if (strict_compare==SPECIAL_TRAINING) strict = "-special";
    ofstream file("regres/lr-train"+strict, fstream::out);
    for(auto d : data)
        For(i,d.size()) file << d[i] << char((i+1==(int)d.size())?'\n':' ');
    file.close();
}
void MachineLinear::load() {
    ifstream file("regres/lr-model", fstream::in);
    int n;
    double x;
    file >> n;
    For(i, n) {
        file >> x;
        coef.push_back(x);
    }
    file >> intercept;
    file.close();
}
void MachineLinearStrict::load() {
    ifstream file("regres/lr-model-strict", fstream::in);
    int n;
    double x;
    file >> n;
    For(i, n) {
        file >> x;
        coef.push_back(x);
    }
    file >> intercept;
    file.close();
}

double MachineLinear::predict(const vector<double>& values) {
    assert(values.size() == coef.size());
    double res = intercept;
    For(i, values.size()) res += values[i]*coef[i];
    //return pow(1./(1.+exp(-res)),4);
    return 1./(1.+exp(-res));
}

size_t calculate_duplication_hash(const Candidate &c, HEvent* event) {
    vector<int> atoms;
    for (auto atom : event->atoms) {
        atoms.push_back(atom.type);
    }

    return boost::hash_value(make_pair(atoms, c.is_inv()));
}

void Machine::add_used_duplication(const Candidate& c, HEvent* event) {
    auto hash = calculate_duplication_hash(c, event);
    used_duplications_now.insert(hash);
}

bool Machine::was_duplication_used(const Candidate& c, HEvent* event) {
    auto hash = calculate_duplication_hash(c, event);
    return used_duplications_prev.find(hash) != used_duplications_prev.end();
}

void Machine::reset_used_duplications() {
    used_duplications_prev.swap(used_duplications_now);
    used_duplications_now.clear();
}
