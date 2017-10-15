#include"sequence.h"

Sequence::Sequence(double age, int length) {
    this->name = "unicorn";
    this->age = age;
    this->first = nullptr;
    if (length > 0) first = new GAtom(length);
}

Sequence::Sequence(Sequence* parent) : Sequence(parent->age, 0) {
    first = (!parent->first)?nullptr:new GAtom(parent->first, parent->first->get_dna());
    ForGAtom(atom, this)
        atom->next = (!atom->next)?nullptr:new GAtom(atom->next, atom->next->get_dna());
}

Sequence::~Sequence() {
    GAtom* atom = first;
    while(atom != nullptr) {
        GAtom* old = atom;
        atom = atom->next;
        delete old;
    }
}

int Sequence::length() {
    int len = 0;
    ForGAtom(atom, this) len+=atom->length();
    return len;
}
int Sequence::com_length() {
    int len = 0;
    ForGAtom(atom, this) len+=atom->com_length();
    return len;
}
int Sequence::atom_count() {
    int len = 0;
    ForGAtom(atom, this) len+=1;
    return len;
}

void Sequence::mutate(double time) {
    age += time;
    ForGAtom(atom, this) atom->mutate(time);
}

void Sequence::split_breakpoints(vector<int> positions) {
    sort(positions.begin(), positions.end());
    GAtom* atom = this->first;
    int pos = 0;
    for(int p : positions) {
        while(atom != nullptr && pos + atom->length() <= p) {
            pos += atom->length();
            atom = atom->next;
        }
        if (pos < p) atom->get_type()->split(atom->is_inverted()?atom->length()-p+pos:p-pos);
    }
}

set<GAtomType*> Sequence::retype_atoms(int length_threshold) {
    set<GAtomType*> res; 
    int name_id = 0;
    ForGAtom(atom, this) {
        if (atom->length() < length_threshold) {
            atom->get_type()->id = 0;
        } else {
            atom->set_name(name.substr(0,1) + to_string(++name_id));
            if (res.count(atom->get_type()) == 0) {
                res.insert(atom->get_type());
                atom->get_type()->id = SIZE(res);
                if (atom->is_inverted()) atom->get_type()->invert();
            }
        }
    }
    return res;
}

void Sequence::write_dna(ostream& os, const string& sep) {
    ForGAtom(atom, this) atom->write_dna(os, "", true);
    os << sep;
}

void Sequence::write_atoms_short(ostream& os, const string& sep) {
    ForGAtom(atom, this) if (atom->get_id()) atom->write_type(os, " ");
    os << sep;
}

void Sequence::write_atoms(ostream& os) {
    int pos = 0;
    ForGAtom(atom, this) {
        if (atom->get_id()) {
            os << name << " " << atom->get_name() << " " << abs(atom->get_id()) << " "
               << sign(atom->get_id()) << " " << pos << " " << pos + atom->com_length() << endl;
        }
        pos += atom->com_length();
    }
}

ostream& operator<<(ostream& os, const Sequence& sequence) {
    os << sequence.age << ": ";
    ForGAtom(atom, &sequence) os << *atom << " ";
    return os << endl;
}

