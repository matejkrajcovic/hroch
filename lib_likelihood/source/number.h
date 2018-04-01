#ifndef NUMBER_H
#define NUMBER_H

#include"constants.h"

class Number{
public:
    double data;
    int zap;
    Number();
    Number(double value);
    Number(double value, bool z);
    void print(FILE* f = stderr);
    void print_double(FILE* f = stdout);
    Number& operator*=(const Number& n);
    Number& operator+=(const Number& n);
    bool zero() const {return data == min_inf;};
    double true_value();
};

bool operator<(const Number& n1, const Number& n2);

Number operator*(const Number& n1, const Number& n2);
Number operator/(const Number& n1, const Number& n2);
Number operator+(const Number& n1, const Number& n2);
Number logNumber(double _data);

#endif
