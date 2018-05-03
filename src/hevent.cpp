#include <cassert>

#include "hevent.h"

using namespace std;

ostream& operator<<(ostream& os, const HEvent& event) {
    os << event.species << " " << event.name << " " << (event.parent?event.parent->name:"root")
       << " " << event.event_time << " " << event.type << " ";

    for(auto a : event.atoms) os << a.type << " ";
    os << "#";
    for(auto a : event.atom_parents) os << " " << a;
    return os << endl;
}

HEvent::HEvent() {
    edge_time = event_time = 0.;
    parent = nullptr;
}

HEvent::HEvent(string species, string name, string type) :
    HEvent::HEvent() {
    this->species = species;
    this->name = name;
    this->type = type;
}

HEvent::HEvent(const string& name, GEvent* event, Sequence* after) :
    HEvent::HEvent() {
    this->species = after->name;
    this->name = name;
    this->edge_time = event->get_time();
    this->event_time = after->age;
    this->type = event->name();
}

HEvent::HEvent(string species, string name, string type, HEvent* same_child) :
    HEvent::HEvent(species, name, type) {
    this->atoms = same_child->atoms;
    same_child->parent = this;
    For(i, same_child->atoms.size()) same_child->atom_parents.push_back(i);
}

HEvent::HEvent(History* history, istringstream& iss) : HEvent::HEvent() {
    string parent_name, word;
    iss >> species >> name >> parent_name >> event_time >> type;
    parent = (parent_name=="root")?nullptr:history->events[parent_name];
    while(iss >> word) {
        if (word == "#") break;
        atoms.push_back(HAtom(stoi(word)));
    }
    while(iss >> word)
        atom_parents.push_back(stoi(word));
    assert(atoms.size() == atom_parents.size());

    if (type == "leaf") atoms = history->leaf_atoms[species];
}

bool HEvent::is_final() {
    set<int> types;
    for(auto& a : atoms) {
        if (types.count(a.atype())) return false;
        types.insert(a.atype());
    }
    type = "root";
    atom_parents.resize(atoms.size());
    For(i, atoms.size()) atom_parents[i] = -1;
    return true;
}

bool HEvent::is_useless() {
    if (type == "dup" || type == "dupi" || type == "del") {
       if (parent != nullptr && atoms.size() == parent->atoms.size())
           return true;
    }
    return false;
}

void HEvent::compute_atoms(HEvent* parent, Sequence* before, Sequence* after) {
    this->parent = parent;
    map<GAtom*, int> position;
    if (before == nullptr) {
        position[nullptr] = -1;
    } else {
        int pos = 0;
        ForGAtom(atom, before) if (atom->get_id()) {
            position[atom] = pos;
            pos++;
        }
    }
    ForGAtom(atom, after) if (atom->get_id()) {
        atoms.push_back(HAtom(atom->get_id()));
        atom_parents.push_back(position[atom->parent]);
    }
}

void HEvent::compute_atom_ids(Sequence* after) {
    int pos = 0;
    ForGAtom(atom, after) if (atom->get_id()) {
        if (atom->get_name().size())
            atoms[pos].add_id(HAtom::str_to_id(atom->get_name()));
        if (parent != nullptr)
            parent->atoms[atom_parents[pos]].add_ids(atoms[pos].get_ids());
        pos++;
    }
    assert(pos == (int) atoms.size());
}

void HEvent::compute_atom_ids(HEvent* after) {
    For(i, after->atoms.size()) {
        atoms[after->atom_parents[i]].add_ids(after->atoms[i].get_ids());
    }
}

int HEvent::compute_is_left() {
    if (type!="dup" && type!="dupi") return -1;
    map<int, int> M;
    for(auto p : atom_parents) M[p]++;
    For(i, atom_parents.size()) if (M[atom_parents[i]] > 1) {
        return (atom_parents[i]==i);
    }
    return -1;
}

void HEvent::compute_diff(const vector<HAtom>& diff) {
    assert(diff.size()==0);
    map<int, int> M;
    for(auto p : atom_parents) M[p]++;
    diff_atoms.clear();
    For(i, atom_parents.size()) if (M[atom_parents[i]] > 1)
        diff_atoms.push_back(atoms[i]);
    sort(diff_atoms.begin(), diff_atoms.end());
}

void HEvent::test_stats(History* h, ostream& os) {
    map<int, vector<HAtom>> M;
    For(i, atom_parents.size()) M[atom_parents[i]].push_back(atoms[i]);
    os << "cheries " << name << " ";
    bool reverse = !compute_is_left();
    for(auto p : M) {
        if (p.second.size() == 2) {
            if (reverse) swap(p.second[0],p.second[1]);
            os << p.second[0].type << "("
               << p.second[0].get_ids() << p.second[1].get_ids()
               << h->cherryness(p.second[0],p.second[1]) << ")  ";
        }
    }
    os << endl;
}

void HEvent::write_detailed(ostream& os) {
    for(const HAtom& atom : atoms) {
        os << atom.type << atom.get_ids() << " ";
    }
    os << endl;
}

void HEvent::clear_parent() {
    if (parent) {
        delete parent;
        parent = nullptr;
        atom_parents.clear();
    }
}
