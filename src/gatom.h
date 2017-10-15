// Representation of atom used when generating sequences

#ifndef ATOM_H
#define ATOM_H

#include"constants.h"

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
    vector<GAtom*> children;
    string dna;
    string name;
    bool inverted;
public:
    GAtom *parent, *next;

    int length() { return dna.size(); }
    int com_length() { return dna.size() - count(dna.begin(), dna.end(), '-'); }
    GAtomType* get_type() { return type; }
    int get_id() { return inverted?-type->id:type->id; }
    const string& get_dna() { return dna; }
    bool is_inverted() { return inverted; }
    string get_name() { return name; } 
    void set_name(const string& str) { name = str; } 

    bool change_parent(GAtom* new_parent); // return if changed
    void unlink_parent() { change_parent(nullptr); }
    void invert();
    void hard_invert();
    void mutate(double time);
    void split(int position, GAtom* first_parent, GAtom* second_parent); 
    GAtom* duplicate();                                      
    GAtom(int length = 0);
    GAtom(GAtom* parent, const string& dna);

    void write_dna(ostream& os = cout, const string& sep = " ", bool compact = false);
    void write_type(ostream& os = cout, const string& sep = " ");
};
ostream& operator<<(ostream& os, GAtom& atom);

#endif
