#ifndef EVENT_H
#define EVENT_H

#include"constants.h"
#include"number.h"
#include"history.h"
#include"atomtree.h"

namespace likelihood {
class Event{
    HistoryLikelihood *history;

    const trint beg = trint(6,0,0), end = trint(6,0,1);
    vector<vvdo> S;
    vector<vvvi> P;
    vector<pair<trint, double> > E;
    void add_edge(int fromk, int fromi, int fromj, int tok, int toi, int toj, double val);
    trint path_previous(trint x);
    void compute_graph();
    void compute_prevgenes(const vtri& path, Event* e, vector<piii>& deletions, vector<int>& data);
    double path_score(const vtri& path);

public:
    string species, name, type;
    double edge_time, event_time;
    int ev_len, ev_dist, sequence_length; //TODO teraz je to v atomoch, zmenit na bazy
    Event* parent;
    vector<int> genes;
    vector<int> gene_parents;
    vector<Node*> atom_nodes;

    Number likelihood();
    void read(istringstream& is);
    void write(ostream& os);
    void writeshort(ostream& os);

    Event(){;};
    Event(HistoryLikelihood* _history);
    Event(HistoryLikelihood* _history, istringstream& is){
        this->history = _history;
        read(is);
    }

    // vyrabanie
    bool final();
    void load_leaf(const pair<string, vector<atom_element> >& atom_list);
    Event* event_from_path(const vtri& path);
    Event* find_prev_event();
    vtri apply_event(const vi&  event);
    vi get_change(); // pre duplikaciu, co sa duplikovalo
};

double best_time(const vector<pair<Node*, Node*> >& nodes);
}
#endif
