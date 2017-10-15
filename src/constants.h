#ifndef CONSTANTS_H
#define CONSTANTS_H

#include<cstdio>
#include<iostream>
#include<sstream>
#include<fstream>
#include<string>
#include<algorithm>
#include<vector>
#include<map>
#include<set>
#include<cmath>
#include<cassert>
#include<queue>
using namespace std;

#define For(i,n) for(int i=0; i<int(n); ++i)
#define ForGAtom(atom, sequence) for(GAtom* atom = (sequence)->first; atom != nullptr; atom = atom->next)
#define SIZE(i) int(i.size())
#define UNUSED(x) (void)(x)
#define BASES 4
#define DATAPATH "data/"
#define epsilon 1e-10

extern int LOWER_RANGE;
extern int UPPER_RANGE;
extern string TEST_CASE;

extern char bases[];
struct trint;
class GAtom;
class GAtomType;
class GEvent;
class Sequence;
class Candidate;
class Dynamics;
class HAtom;
class HEvent;
class History;
class GHistory;
class CherryForest;
class Machine;
class ScoringData;

typedef pair<int, int> pii;
typedef pair<pii, int> piii;
typedef pair<pii, double> piid;

typedef vector<pii> vpi;
typedef vector<vpi> vvpi;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vvi> vvvi;
typedef vector<double> vdo;
typedef vector<vdo> vvdo;
typedef vector<trint> vtri;

#include"random.h"
#include"gatom.h"
#include"gevent.h"
#include"hatom.h"
#include"hevent.h"
#include"cherry.h"
#include"sequence.h"
#include"history.h"
#include"ghistory.h"
#include"dynamics.h"
#include"score.h"
#include"machine.h"

struct trint {
    int a,b,c;
    trint() {a = b = c = 0;}
    trint(const int& a, const int& b, const int& c) {
        this->a = a;
        this->b = b;
        this->c = c;
    }
};
inline bool operator==(const trint& t1, const trint& t2) {
    return t1.a==t2.a && t1.b==t2.b && t1.c==t2.c;
}
inline bool operator!=(const trint& t1, const trint& t2) {
    return t1.a!=t2.a || t1.b!=t2.b || t1.c!=t2.c;
}
inline bool operator<(const trint& t1, const trint& t2) {
    if (t1.a!=t2.a) return t1.a<t2.a;
    if (t1.b!=t2.b) return t1.b<t2.b;
    return t1.c<t2.c;
}

template<typename T>
ostream& operator<<(ostream& os, const vector<T> v) {
    os << "[";
    For(i, SIZE(v)) (i?(os << ", "):os) << v[i]; 
    return os << ']';
}
template <typename T> 
int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

extern int fail_on_error;
extern int error_happened;
extern int debuging;
extern int stats;
extern int do_cheeryness;
#define SPECIAL_TRAINING 47
extern int strict_compare;
extern char bases[BASES];
extern int base_id[256];
extern char base_inv[256];

inline void mustbe(const bool& True, const string& message="") {
    if (True!=true) {
        cerr << "Error: " << message << endl; 
        exit(1);
    }
}

void setup_constants();
set<string> parse_arguments(int argc, char **argv);

#endif
