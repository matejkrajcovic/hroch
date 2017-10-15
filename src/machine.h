#ifndef MACHINE_H
#define MACHINE_H

#include"constants.h"

class Machine {
public:
    virtual void train_data(const vector<double>& values, double result) = 0;
    virtual void save() {;}
    virtual void load() {;}
    virtual double predict(const vector<double>& values) = 0;
    virtual ~Machine() {};
};

class MachineOne: public Machine {
    virtual void train_data(const vector<double>& values, double result) {UNUSED(values); UNUSED(result);}
    virtual double predict(const vector<double>& values) {UNUSED(values); return 1.0;}
};

class MachineBachelor: public MachineOne {
public:
    double predict(const vector<double>& values);
};

class MachineLinear: public Machine {
    vector<vector<double>> data;
protected:
    vector<double> coef;
    double intercept;
public:
    void train_data(const vector<double>& values, double result);
    void save();
    void load();
    double predict(const vector<double>& values);
    
};

class MachineLinearStrict: public MachineLinear {
public:
    void load();
};

#endif
