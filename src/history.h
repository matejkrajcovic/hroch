#ifndef HISTORY_H
#define HISTORY_H

#include"constants.h"

// event similarity
#define DIFFERENT 0
#define SAME_ANY 1
#define SAME_LAST 2

// scoring strategy
#define NO_STRATEGY 0
#define CHERRY_NO 1
#define CHERRY_TREE 2
#define CHERRY_LEN 3
#define SCORE_CL 4
#define SCORE_BAC_NC 5
#define SCORE_BAC 6
#define SCORE_LR 7
#define SCORE_LRS 8
#define KNOW_HOW 9

#define DEFAULT_CHERRY CHERRY_TREE

// evaluation
#define EVAL_INF -1
#define EVAL_LAZY -2

class History {
    int strategy, cherry_mode;
    Machine* machine;
    void init_zero();
    vector<string> leaf_species;
    map<string, string> leaf_dna_sequence;
    map<int, string> leaf_atom_dna;
    CherryForest* cherry_forest;
public:
    History* original;
    map<string, double> stats;
    map<string, vector<HAtom>> leaf_atoms;
    map<string, HEvent*> events; // name is key
    map<string, HEvent*> leaf_events; // species is key
    map<int, int> type_dna_length;

    void set_strategy(int strategy, Machine* machine = nullptr);

    void proc_learn();
    void proc_reconstruct(int number = EVAL_INF);
    void real_reconstruct();
    void proc_test_candi(int strategy, string mark);
    void proc_test_score(int strategy, string mark);
    void rec_parent(HEvent* event);
    void rec_compute_parent(const Candidate& C, HEvent* event);
    HEvent* rec_see_event(const Candidate& C, HEvent* event);
    void rec_setup_scoring_data(const Candidate& C, HEvent* event, ScoringData* sd);
    void rec_merge_candidate(const Candidate& C, HEvent* event);
    set<Candidate> rec_candidates(HEvent* event);
    double rec_score(const Candidate& c, HEvent* event);

    HEvent* nth_from_end(int n);
    HEvent* resolve_deletion(HEvent* deletion);
    int is_original(HEvent* event, bool strict = false);
    int is_correct(bool weak = false);
    double cherryness(const HAtom& a, const HAtom& b, int mode = 0);
    void merge(const HAtom& a, const HAtom& b);

    void read_final_sequence(string species, istream& is);
    void read_atoms(istream& is);
    void read_atoms_align(const string& basepath);
    void read_events(istream& is);
    void read_cherryness(const string& basepath);

    string gen_event_name();
    void save(string name);
    void write_events(ostream& os);
    void write_stats(ostream& os);
    double get_time();
    bool same_as(History *h);

    History(string atoms_file, string trees_dir, int strategy);
    History(string basename, string id = "");
    History(History* original);
    ~History();
};


#endif
