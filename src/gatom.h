// Representation of atom used when generating sequences

#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>
#include <iostream>
#include <ostream>
#include <algorithm>

class GAtom;

class GAtomType {
    static int type_cnt;
public:
    static void reset_ids() {type_cnt = 0;}
    int id;
    GAtom* root;

    void split(int position);
    void invert();
    GAtomType(GAtom* root);
};

class GAtom {
    GAtomType* type;
    std::vector<GAtom*> children;
    std::string dna;
    std::string name;
    bool inverted;
public:
    GAtom *parent, *next;

    int length() { return dna.size(); }
    int com_length() { return dna.size() - std::count(dna.begin(), dna.end(), '-'); }
    GAtomType* get_type() { return type; }
    int get_id() { return inverted?-type->id:type->id; }
    const std::string& get_dna() { return dna; }
    bool is_inverted() { return inverted; }
    std::string get_name() { return name; }
    void set_name(const std::string& str) { name = str; }

    bool change_parent(GAtom* new_parent); // return if changed
    void unlink_parent() { change_parent(nullptr); }
    void invert();
    void hard_invert();
    void mutate(double time);
    void split(int position, GAtom* first_parent, GAtom* second_parent);
    GAtom* duplicate();
    GAtom(int length = 0);
    GAtom(GAtom* parent, const std::string& dna);

    void write_dna(std::ostream& os = std::cout, const std::string& sep = " ", bool compact = false);
    void write_type(std::ostream& os = std::cout, const std::string& sep = " ");
};

std::ostream& operator<<(std::ostream& os, GAtom& atom);

#endif
