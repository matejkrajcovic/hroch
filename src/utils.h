#ifndef UTILS_H
#define UTILS_H

#include <ostream>
#include <vector>
#include <algorithm>

#define For(i,n) for(int i=0; i<int(n); ++i)
#define ForGAtom(atom, sequence) for(GAtom* atom = (sequence)->first; atom != nullptr; atom = atom->next)
#define UNUSED(x) (void)(x)

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> v) {
    if (!v.size()) {
        os << "[]";
        return os;
    }
    os << "[";
    std::for_each(v.begin(), --v.end(), [&](T elem){
        os << elem << ", ";
    });
    os << v.back() << "]";
    return os;
}

template <typename T>
int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif
