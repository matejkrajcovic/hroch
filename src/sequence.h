// Representation of DNA sequence used for generated data

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <ostream>

#include "gatom.h"

class Sequence {
public:
    GAtom* first;
    double age;
    std::string name;

    Sequence(double age = 0.0, int length = 0);
    Sequence(Sequence* parent);
    ~Sequence();
    int length();
    int com_length();
    int atom_count();
    void mutate(double time);
    void split_breakpoints(std::vector<int> positions);

    std::set<GAtomType*> retype_atoms(int length_threshold = 0);
    void write_dna(std::ostream& os = std::cout, const std::string& sep = "\n");
    void write_atoms_short(std::ostream& os = std::cout, const std::string& sep = "\n");
    void write_atoms(std::ostream& os = std::cout);
};

std::ostream& operator<<(std::ostream& os, const Sequence& sequence);

#endif
