#include <string>
#include <fstream>
#include <cmath>
#include <cassert>
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

void Machine::add_used_duplication(HEvent* event) {
    HEvent* parent = event;
    while (parent->parent) {
        parent = parent->parent;
    }
    auto slices = get_changed_slices_of_event(event, parent);
    for (auto slice : slices) {
        used_duplications_now.insert(slice);
        used_duplications_now.insert(get_inverse_slice(slice));
    }
}

bool Machine::was_duplication_used(HEvent* event) {
    HEvent* parent = event;
    while (parent->parent) {
        parent = parent->parent;
    }
    auto slices = get_changed_slices_of_event(event, parent);
    bool all_used = true;
    for (auto slice : slices) {
        if (used_duplications_prev.find(slice) == used_duplications_prev.end()) {
            all_used = false;
            break;
        }
    }

    return all_used;
}

void Machine::reset_used_duplications() {
    used_duplications_prev.swap(used_duplications_now);
    used_duplications_now.clear();
}
