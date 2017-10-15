// Representation of DNA sequence used for generated data

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include"constants.h"

class Sequence {
public:
    GAtom* first;
    double age;
    string name;
   
    Sequence(double age = 0.0, int length = 0);
    Sequence(Sequence* parent);
    ~Sequence();
    int length();
    int com_length();
    int atom_count();
    void mutate(double time);
    void split_breakpoints(vector<int> positions);

    set<GAtomType*> retype_atoms(int length_threshold = 0);
    void write_dna(ostream& os = cout, const string& sep = "\n");
    void write_atoms_short(ostream& os = cout, const string& sep = "\n");
    void write_atoms(ostream& os = cout);
};

ostream& operator<<(ostream& os, const Sequence& sequence);

#endif
