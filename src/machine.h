#ifndef MACHINE_H
#define MACHINE_H

#include <vector>
#include <set>

#include "candidate.h"
#include "utils.h"

class HEvent;

class Machine {
public:
    std::set<std::vector<int>> used_duplications_now;
    std::set<std::vector<int>> used_duplications_prev;
    virtual void train_data(const std::vector<double>& values, double result) = 0;
    virtual void save() {;}
    virtual void load() {;}
    virtual double predict(const std::vector<double>& values) = 0;
    virtual ~Machine() {};
    void add_used_duplication(HEvent* event);
    bool was_duplication_used(HEvent* event);
    void reset_used_duplications();
};

class MachineOne: public Machine {
    virtual void train_data(const std::vector<double>& values, double result) {UNUSED(values); UNUSED(result);}
    virtual double predict(const std::vector<double>& values) {UNUSED(values); return 1.0;}
};

class MachineBachelor: public MachineOne {
public:
    double predict(const std::vector<double>& values);
};

class MachineLinear: public Machine {
    std::vector<std::vector<double>> data;
protected:
    std::vector<double> coef;
    double intercept;
public:
    void train_data(const std::vector<double>& values, double result);
    void save();
    void load();
    double predict(const std::vector<double>& values);

};

class MachineLinearStrict: public MachineLinear {
public:
    void load();
};

#endif
