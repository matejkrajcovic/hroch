#include"number.h"

Number::Number(){
    data = min_inf;
    zap = 0;
}

Number::Number(double value){
    data = log(abs(value));
    if(value<0) zap=1;
    else zap = 0;
}

Number logNumber(double _data){
    Number n;
    n.data = _data;
    n.zap = 0;
    return n;
}

Number::Number(double value, bool z) {
    if(z) data = min_inf;
    else data = log(value);
    if(value < 0) zap=1;
    else zap=0;
}

void Number::print(FILE* f){
    fprintf(f,"%c %.12lf\n",(this->zap == 0)?'+':'-' ,this->data);
}

void Number::print_double(FILE* f){
    fprintf(f,"%.20lf\n", exp(this->data));
}

Number& Number::operator*=(const Number& n){
    if(n.zero()) {data = min_inf; zap=0;}
    else if(!zero()) {
        data += n.data;
        zap = zap * n.zap;
    }
    return *this;
}

Number& Number::operator+=(const Number& n) {
    if(n.zero()) return *this;
    if(zero()) {data = n.data; zap = n.zap; return *this;}
    if(zap == n.zap) {
        if(data < n.data) data = n.data + log(1 + exp(data - n.data));
        else data = data + log(1 + exp(n.data - data));
    }
    else {
        if(data < n.data) {
            zap = n.zap;
            data = n.data + log(1 - exp(data - n.data));
        }
        else {
            data = data + log(1 - exp(n.data - data));
        }
    }
    return *this;
}

bool operator<(const Number& n1, const Number& n2){
    if(n1.zap != n2.zap) {
        if(n1.zap==1) return true;
        return false;
    }
    if(n1.zap == 1) return (n1.data > n2.data);
    return (n1.data < n2.data);
}

Number operator*(const Number& n1, const Number& n2){
    if(n1.zero()) return Number();
    if(n2.zero()) return Number();
    Number res = Number();
    res.data = n1.data + n2.data;
    res.zap = n1.zap * n2.zap;
    return res;
}

Number operator/(const Number& n1, const Number& n2){
    if(n1.zero()) return Number();
    Number res = Number();
    res.zap = n1.zap * n2.zap;
    res.data = n1.data - n2.data;
    return res;
}

Number operator+(const Number& n1, const Number& n2) {
    if(n1.zero()) return n2;
    if(n2.zero()) return n1;
    if(n1.zap == n2.zap) {
        Number res = Number();
        res.zap = n1.zap;
        if(n1.data < n2.data) res.data = n2.data + log(1+exp(n1.data - n2.data));
        else res.data = n1.data + log(1+exp(n2.data - n1.data));
        return res;
    }
    else {
        Number res = Number();
        if(n1.data < n2.data) {
            res.zap = n2.zap;
            res.data = n2.data + log(1 - exp(n1.data - n2.data));
        }
        else {
            res.zap = n1.zap;
            res.data = n1.data + log(1 - exp(n2.data - n1.data));
        }
        return res;
    }
}

double Number::true_value() {
    double res = exp(data);
    if(zap==1) res*=-1;
    return res;
}
