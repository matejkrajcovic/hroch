// Representation of events during reconstruction

#ifndef HEVENT_H
#define HEVENT_H

#include <vector>
#include <ostream>
#include <sstream>

#include "sequence.h"
#include "gevent.h"
#include "history.h"
#include "hatom.h"

class History;

class HEvent {
public:
    std::vector<HAtom> diff_atoms;

    std::string species, name, type;
    double edge_time, event_time;
    HEvent* parent;

    std::vector<HAtom> atoms;
    std::vector<int> atom_parents;

    HEvent();
    HEvent(std::string species, std::string name, std::string type);
    HEvent(std::string species, std::string name, std::string type, HEvent* same_child);
    HEvent(History* history, std::istringstream& iss);
    HEvent(const std::string& name, GEvent* event, Sequence* after);

    bool is_final();
    bool is_useless();

    // call when you know the atoms
    void compute_atoms(HEvent* parent, Sequence* before, Sequence* after);
    // call from leaves to root
    void compute_atom_ids(Sequence* after);
    void compute_atom_ids(HEvent* after);
    // ...
    void compute_diff(const std::vector<HAtom>& diff = {});
    int compute_is_left();
    void test_stats(History* h, std::ostream& os);
    void write_detailed(std::ostream& os);

    friend bool operator==(const HEvent& e1, const HEvent& e2) {
        if (e1.type != e2.type) return false;
        if (e1.type == "root" || e1.type == "leaf") return true;
        return e1.diff_atoms == e2.diff_atoms;
    }
};

std::ostream& operator<<(std::ostream& os, const HEvent& event);

#endif
