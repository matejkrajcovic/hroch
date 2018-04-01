#ifndef HISTORY_H
#define HISTORY_H

#include"constants.h"
#include"number.h"
#include"event.h"
#include"atomtree.h"

// Tato trieda bude zdruzovat informacie pre historiu a jej pomocne veci
// bude obsahovat:
//      Pravdepodobnostne parametre
//      Strom eventov
//      Zarovnania atomov
//      Zoznam atomov v sekvenciach druhov
//      Atomove stromceky


class History{
    // znaci si, co uz pozna
    bool kn_parameters, kn_atom_lists, kn_events, kn_trees;
    int kn_alignments;
public:
    map<string, Event*> events;
    map<string, string> atom_string;                                            // dna sekvencia atomu
    map<string, int> atom_id;                                                   // typ atomu
    map<int, vector<string> > group_atoms;                                      // zoznam atomov pre typ
    map<int, AtomTree*> group_trees;                                             // atomovy strom pre typ
    map<string, double> parameters;
    map<string, vector<atom_element> > atom_list;                               // postupnost atomov pre druh
    set<int> existing_atoms;                                                    // mnozina typov
     
    void read_parameters(istream& is);
    void read_atom_lists(istream& is);
    void read_atom_alignments(istream& is, int group_id);
    void read_events(istream& is);
    void write_events(ostream& os);

    void compute_events();
    void clear_events();

    // vyrobí strom eventov a atómové stromčeky
    void build_trees();

    Number likelihood_all();
    History();

    void trust_me(){
        kn_events = true;
    }
    
    //upravovanie dlzky hran historie
    void edit_edge_lengths();
    void find_better_time(Event* edit, double max_time);
    //ternarne upravovanie dlzky
    void find_better_time_ternary(Event* edit, double max_time);
    Number ternary_full_likelihood(Event* edit, double max_time, double time);
    
    // tu si pamatame zaujimave veci zo starej historie
    map<vi, double> old_events, old_seqences; 
    void add_to_old(Event* e, double koef);

    // zapnute skorovanie
    int scoring_length, scoring_similar, scoring_cherryness, 
        scoring_boundary, scoring_oldhistory, felsenstein_on,
        scoring_divisor;

    // ine
    double branch_time = 0.04;
};

#endif
