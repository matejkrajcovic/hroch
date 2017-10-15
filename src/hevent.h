// Representation of events during reconstruction

#ifndef HEVENT_H
#define HEVENT_H

#include<algorithm>
#include<vector>
using namespace std;

class Candidate {
public:
    int b1,e1,b2,e2;
    inline bool is_inv() const { return e2 < b2; }
    inline bool is1_left() const { return e1 < e2; }
    inline bool is1_right() const { return e2 < e1; }
    inline bool is_valid() const { return e1<min(b2+1,e2) || max(b2-1,e2)<b1+1; }
    vector<int> directions;
    Candidate(int b1, int e1, int b2, int e2);
    Candidate(int b1, int e1, int b2, int e2, const vector<int>& directions);

    void swap_dir();
    
    friend bool operator==(const Candidate& c1, const Candidate& c2) {
        return (c1.b1 == c2.b1) && 
               (c1.e1 == c2.e1) &&
               (c1.b2 == c2.b2) &&
               (c1.e2 == c2.e2) &&
               (c1.directions == c2.directions);
    }
    friend bool operator<(const Candidate& c1, const Candidate& c2) {
        if (c1.b1 != c2.b1) return c1.b1 < c2.b1;
        if (c1.e1 != c2.e1) return c1.e1 < c2.e1;
        if (c1.b2 != c2.b2) return c1.b2 < c2.b2;
        if (c1.e2 != c2.e2) return c1.e2 < c2.e2;
        return c1.directions < c2.directions;
    }
};

#include"constants.h"

ostream& operator<<(ostream& os, const Candidate& c);

class HEvent {
public:
    vector<HAtom> diff_atoms;
    
    string species, name, type;
    double edge_time, event_time;
    HEvent* parent;
  
    vector<HAtom> atoms;
    vector<int> atom_parents;

    HEvent();
    HEvent(string species, string name, string type);
    HEvent(string species, string name, string type, HEvent* same_child);
    HEvent(History* history, istringstream& iss);
    HEvent(const string& name, GEvent* event, Sequence* after);

    bool is_final();
    bool is_useless();

    // call when you know the atoms
    void compute_atoms(HEvent* parent, Sequence* before, Sequence* after);
    // call from leaves to root
    void compute_atom_ids(Sequence* after);
    void compute_atom_ids(HEvent* after);
    // ...
    void compute_diff(const vector<HAtom>& diff = {});
    int compute_is_left();
    void test_stats(History* h, ostream& os);
    void write_detailed(ostream& os);

    friend bool operator==(const HEvent& e1, const HEvent& e2) {
        if (e1.type != e2.type) return false;
        if (e1.type == "root" || e1.type == "leaf") return true;
        return e1.diff_atoms == e2.diff_atoms;
    }
};

ostream& operator<<(ostream& os, const HEvent& event);

#endif
