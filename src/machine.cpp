#include"machine.h"

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
        For(i,SIZE(d)) file << d[i] << char((i+1==SIZE(d))?'\n':' ');
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
    assert(SIZE(values) == SIZE(coef));
    double res = intercept;
    For(i, SIZE(values)) res += values[i]*coef[i];
    //return pow(1./(1.+exp(-res)),4); 
    return 1./(1.+exp(-res)); 
}
