#ifndef CONSTANTS_H
#define CONSTANTS_H

#include<cstdio>
#include<iostream>
#include<sstream>
#include<string>
#include<algorithm>
#include<vector>
#include<map>
#include<set>
#include<cmath>
#include<cassert>
using namespace std;

#include"random.h"

#define For(i,n) for(typeof(n) i=0; i<(n); ++i)
#define SIZE(i) int(i.size())
const double min_inf = log(0);

extern char bases[];
class AtomTree;
class Event;
class History;

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
typedef pair<string,int> atom_element;

struct trint{
    int a,b,c;
    trint(){a = b = c = 0;}
    trint(const int& a, const int& b, const int& c){
        this->a = a;
        this->b = b;
        this->c = c;
    }
};
inline bool operator==(const trint& t1, const trint& t2){
    return t1.a==t2.a && t1.b==t2.b && t1.c==t2.c;
}
inline bool operator!=(const trint& t1, const trint& t2){
    return t1.a!=t2.a || t1.b!=t2.b || t1.c!=t2.c;
}
inline bool operator<(const trint& t1, const trint& t2){
    if (t1.a!=t2.a) return t1.a<t2.a;
    if (t1.b!=t2.b) return t1.b<t2.b;
    return t1.c<t2.c;
}

typedef vector<trint> vtri;

extern int debuging;

#define FEL_ALPHA (1.0/3.0)

inline void mustbe(const bool& True, const string& message=""){
    if (True!=true){
        cerr << "Error: " << message << endl; 
        exit(1);
    }
}

#endif
