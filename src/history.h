#ifndef HISTORY_H
#define HISTORY_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <istream>
#include <ostream>

#include "machine.h"
#include "cherry.h"
#include "hatom.h"
#include "hevent.h"
#include "candidate.h"
#include "score.h"

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

class CherryForest;
class ScoringData;

class History {
    int strategy, cherry_mode;
    Machine* machine;
    void init_zero();
    std::vector<std::string> leaf_species;
    std::map<std::string, std::string> leaf_dna_sequence;
    std::map<int, std::string> leaf_atom_dna;
    CherryForest* cherry_forest;
public:
    History* original;
    std::map<std::string, double> stats;
    std::map<std::string, std::vector<HAtom>> leaf_atoms;
    std::map<std::string, HEvent*> events; // name is key
    std::map<std::string, HEvent*> leaf_events; // species is key
    std::map<int, int> type_dna_length;

    void set_strategy(int strategy, Machine* machine = nullptr);

    void proc_learn();
    void proc_reconstruct(int number = EVAL_INF);
    bool real_reconstruct();
    void proc_test_candi(int strategy, std::string mark);
    void proc_test_score(int strategy, std::string mark);
    void rec_parent(HEvent* event);
    void rec_compute_parent(const Candidate& C, HEvent* event);
    HEvent* rec_see_event(const Candidate& C, HEvent* event);
    void rec_setup_scoring_data(const Candidate& C, HEvent* event, ScoringData* sd);
    void rec_merge_candidate(const Candidate& C, HEvent* event);
    std::set<Candidate> rec_candidates(HEvent* event);
    double rec_score(const Candidate& c, HEvent* event);

    int get_history_score_num_events();
    double get_history_score_likelihood(std::string atoms_filename, std::string align_dir);

    std::set<std::vector<int>> get_changed_slices(bool dels_only_in_dups = true);

    HEvent* nth_from_end(int n);
    HEvent* resolve_deletion(HEvent* deletion);
    int is_original(HEvent* event, bool strict = false);
    int is_correct(bool weak = false);
    double cherryness(const HAtom& a, const HAtom& b, int mode = 0);
    void merge(const HAtom& a, const HAtom& b);

    void read_final_sequence(std::string species, std::istream& is);
    void read_atoms(std::istream& is);
    void read_atoms_align(const std::string& basepath);
    void read_events(std::istream& is);
    void read_cherryness(const std::string& basepath);

    std::string gen_event_name();
    void save(std::string name);
    void write_events(std::ostream& os);
    void write_stats(std::ostream& os);
    double get_time();
    bool same_as(History *h);

    History(std::string atoms_file, std::string trees_dir, int strategy);
    History(std::string basename, std::string id = "");
    History(History* original);
    ~History();
};

double calculate_jaccard_index(std::set<std::vector<int>> A, std::set<std::vector<int>> B);

#endif
