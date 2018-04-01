#include "function.h"

Function::Function() {
    data.resize(1, Number());
    min_pow = 0;
}

Function::Function(double value) {
    data.resize(1, Number(value));
    min_pow = 0;
}

Function::Function(double value, bool z) {
    min_pow = 0;
    if(z) data.resize(1,Number(value,true));
    else data.resize(1,Number(value));
}

Function::Function(vector<Number> value, int min_pow_) {
    data = value;
    min_pow = min_pow_;
}

Function::Function(double value1, double value2, int min_pow_) {
    data.resize(2); data[0]=Number(value1); data[1]=Number(value2);
    min_pow = min_pow_;
}

Function Function::derivate() {
    int size = this->data.size();
    vector<Number> res; res.resize(size);
    For(i, this->data.size())
        res[i] = this->data[i] * Number(4 * FEL_ALPHA * (this->min_pow + i));
    return Function(res, this->min_pow);
}

Number Function::evaluate(double time) {
    Number res = Number();
    For(i,this->data.size())
        res += this->data[i] * logNumber(4*FEL_ALPHA * (this->min_pow + i) * time);
    return res;
}

Function operator+(const Function& f1, const Function& f2) {
    int new_min_pow = min(f1.min_pow, f2.min_pow);
    int max_pow = max(f1.min_pow + f1.data.size() - 1, f2.min_pow + f2.data.size() - 1);
    int size = max_pow - new_min_pow + 1;
    vector<Number> res; res.resize(size, Number());
    For(i,f1.data.size()) res[f1.min_pow + i - new_min_pow] += f1.data[i];
    For(i,f2.data.size()) res[f2.min_pow + i - new_min_pow] += f2.data[i];
    return Function(res, new_min_pow);
}

Function operator*(const Function& f1, const Function& f2) {
    int new_min_pow = f1.min_pow + f2.min_pow;
    int max_pow = f1.min_pow + f1.data.size() + f2.min_pow + f2.data.size() -2;
    int size = max_pow - new_min_pow + 1;
    vector<Number> res; res.resize(size, Number());
    For(i,f1.data.size())
        For(j,f2.data.size()) {
            res[f1.min_pow + f2.min_pow + i + j - new_min_pow] += f1.data[i] * f2.data[j];
        }
    return Function(res, new_min_pow);
}

/*void Function::optimize() {
    int z=0;
    while(abs(data[z])<0.000000000001) z++;
    int k=data.size()-1;
    while(abs(data[k])<0.000000000001) k--;
    k++;
    min_pow += z;
    vdo data1;
    for(int i=z; i<k; i++) data1.push_back(data[i]);
    data = data1;
}*/

void Function::print() {
    For(i,data.size()) fprintf(stderr,"%10d ",min_pow+i); fprintf(stderr,"\n");
    For(i,data.size()) data[i].print(); fprintf(stderr,"\n");
}
