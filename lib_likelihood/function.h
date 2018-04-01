#ifndef FUNCTION_H
#define FUNCTION_H

#include "constants.h"
#include "number.h"

namespace likelihood {
class Function {
public:
    vector<Number> data;
    int min_pow;
    Function();
    Function(double value);
    Function(double value, bool z);
    Function(double value1, double value2, int min_pow_);
    Function(vector<Number> value, int min_pow_);
    Function derivate();
    Number evaluate(double time);
    void print();
    //void optimize();
};

Function operator+(const Function& f1, const Function& f2);
Function operator*(const Function& f1, const Function& f2);
}
#endif
