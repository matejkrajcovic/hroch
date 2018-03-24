#include <fstream>
#include <cassert>
#include <cmath>

#include "history.h"

using namespace std;

void open_check(ifstream &f, string filename) {
    f.open(filename, fstream::in);
    if (!f.is_open()) {
        cerr << "Error opening file " << filename << endl;
        exit(1);
    }
}

void History::init_zero() {
    machine = nullptr;
    cherry_mode = DEFAULT_CHERRY;
    strategy = NO_STRATEGY;
}

History::History(string atoms_file, string trees_dir, int strategy) {
    assert(strategy != 0);
    init_zero();
    ifstream f;
    open_check(f, atoms_file);
    read_atoms(f);
    f.close();

    cherry_forest = nullptr;
    if (do_cherryness) {
        read_cherryness(trees_dir + "/");
    }

    for(auto la : leaf_atoms) {
        HEvent* e = new HEvent(la.first, gen_event_name(), "leaf");
        e->atoms = la.second;
        events[e->name] = e;
        leaf_events[e->species] = e;
    }
}

History::History(string basename, string id) {
    init_zero();
    if (id.size()) basename += "-" + id;
    ifstream f;
    cout << "Loading " << basename << endl;

    for(string species : {"unicorn"}) {
        open_check(f, basename+"-"+species+".dna");
        read_final_sequence(species, f);
        f.close();
    }

    open_check(f, basename+".atoms");
    read_atoms(f);
    f.close();

    // We actualt dont need the alignments
    /* read_atoms_align(basename+"/"); */

    open_check(f, basename+".nhistory");
    read_events(f);
    f.close();
    cherry_forest = nullptr;
    if (do_cherryness) {
        for(string species : {"unicorn"}) {
            read_cherryness(datapath + "dupstemp_data-generated-" + id + "-" + species + "-dna/");
        }
    }
    cout << "       " << basename << " loaded" << endl;
}

string History::gen_event_name() {
    return "e-"+to_string(events.size());
}

History::History(History* original) {
    init_zero();
    this->original = original;
    this->leaf_species = original->leaf_species;
    this->leaf_atoms = original->leaf_atoms;
    this->type_dna_length = original->type_dna_length;
    for(auto la : leaf_atoms) {
        HEvent* e = new HEvent(la.first, gen_event_name(), "leaf");
        e->atoms = la.second;
        events[e->name] = e;
        leaf_events[e->species] = e;
    }
    if (do_cherryness)
        this->cherry_forest = new CherryForest(original->cherry_forest);
}

History::~History() {
    for(auto e : events) delete e.second;
    events.clear();
    leaf_events.clear();
    leaf_species.clear();
    leaf_atoms.clear();
    leaf_atom_dna.clear();
    if (cherry_forest != nullptr) delete cherry_forest;
}


void History::read_final_sequence(string species, istream& is) {
    string word;
    is >> word >> leaf_dna_sequence[species];
    assert(">"+species == word);
}

void History::read_atoms(istream& is) {
    string species, name;
    int type, inverz, from, to;
    while(is >> species >> name >> type >> inverz >> from >> to) {
        type *= inverz;
        leaf_atoms[species].push_back(HAtom(type, HAtom::str_to_id(name)));
        type_dna_length[type] = max(type_dna_length[type], to-from);
        type_dna_length[-type] = max(type_dna_length[-type], to-from);
    }
    for(auto la : leaf_atoms) leaf_species.push_back(la.first);
}

void History::read_atoms_align(const string& basepath) {
    map<int, ifstream> files;
    for(auto atoms : leaf_atoms) for(HAtom atom : atoms.second)
        if (files.count(atom.atype()) == 0)
            open_check(files[atom.atype()],
                       basepath+to_string(atom.atype())+".aln");
    for(auto& f : files) {
        string word, buffer, name;
        while(f.second >> word) {
            if (word[0] == '>') {
                if (name.size()) leaf_atom_dna[HAtom::str_to_id(name)] = buffer;
                buffer.clear();
                name = word.substr(1);
            } else buffer += word;
        }
        if (name.size()) leaf_atom_dna[HAtom::str_to_id(name)] = buffer;
        f.second.close();
    }
}

void History::read_events(istream& is) {
    string line;
    HEvent* e = nullptr;
    while(getline(is, line)){
        istringstream iss(line);
        e = new HEvent(this, iss);
        events[e->name] = e;
        if (e->type == "leaf")
            leaf_events[e->species] = e;
    }
    e->atoms = leaf_atoms[e->species];
    while(e->parent != nullptr) {
        e->parent->compute_atom_ids(e);
        e->compute_diff();
        e = e->parent;
    }
    // remove useless events
    auto temp_events = events;
    for(auto ev : temp_events) {
        while ((!ev.second->is_useless()) && ev.second->parent != nullptr &&
            ev.second->parent->is_useless()) {
            ev.second->parent = ev.second->parent->parent;
        }
    }
    for(auto ev : temp_events) {
        if (ev.second->is_useless()) {
            events.erase(ev.first);
            delete ev.second;
        }
    }
}

void History::read_cherryness(const string& basepath) {
    cherry_forest = new CherryForest();
    set<int> atom_types;
    for(auto la : leaf_atoms) for (auto atom : la.second)
        atom_types.insert(atom.atype());
    for(const auto& at : atom_types) {
        cherry_forest->read_atom(this, at, basepath + to_string(at) + ".nex.trprobs");
    }
}

bool compare(HEvent* e1, HEvent* e2) {
    return e1->event_time < e2->event_time;
}

void History::write_events(ostream& os) {
    vector<HEvent*> just_events;
    for(auto e : events) just_events.push_back(e.second);
    sort(just_events.begin(), just_events.end(), compare);
    for(auto e : just_events) os << *e;
}

void History::save(string name) {
    ofstream file("outputs/"+name);
    write_events(file);
    file.close();
}

void History::write_stats(ostream& os) {
    assert(leaf_species.size()==1);
    string species = leaf_species[0];
    os << species << " " << leaf_atoms[species].size() << " atoms ";
    os << original->events.size() << " events; time: " << original->get_time();
    os << " (dc" << do_cherryness << ")" << endl;
    original->nth_from_end(1)->test_stats(this, os);

    for(auto sp : stats) {
        os << sp.first << " " << sp.second << endl;
    }
}

double History::get_time() {
    double min_time = 0.0, max_time = 0.0;
    for(auto e : events) {
        min_time = min(min_time, e.second->event_time);
        max_time = max(max_time, e.second->event_time);
    }
    return max_time - min_time;
}

int History::is_original(HEvent* event, bool strict) {
    assert(event->atom_parents.size());
    event->compute_diff();
    for(auto ev : original->events)
        if (*(ev.second) == *(event)) {
            if (strict && ev.second->compute_is_left() !=
                event->compute_is_left()) continue;

            return (original->nth_from_end(1)==ev.second)?SAME_LAST:SAME_ANY;
        }
    return DIFFERENT;
}

bool History::same_as(History *h) {
    if (events.size() != h->events.size()) return false;
    for(auto ev : events) ev.second->compute_diff();
    for(auto ev : h->events) ev.second->compute_diff();
    for(auto ev : events) {
        bool found = false;
        for(auto eu : h->events) if(*(ev.second) == *(eu.second)) found = true;
        if (! found) return false;
    }
    return true;
}

int History::is_correct(bool weak) {
    if (!weak && events.size() != original->events.size())
        return false;
    for(auto ev : events)
        if(!is_original(ev.second)) {
            return false;
        }
    return true;
}

HEvent* History::nth_from_end(int n) {
    assert(leaf_events.size()==1);
    HEvent* last = leaf_events.begin()->second;
    For(i, n) last = last->parent;
    return last;
}

HEvent* History::resolve_deletion(HEvent* deletion) {
    assert(deletion->parent != nullptr);
    vector<int> ptc;
    set<int> delp;
    for(int p : deletion->atom_parents) delp.insert(p);
    HEvent* cdup = deletion->parent;
    For(i, cdup->atoms.size()) if (!delp.count(i)) ptc.push_back(i);

    while(cdup != nullptr) {
        if (cdup->type == "dup" || cdup->type == "dupi") {
            bool inside = false;
            bool border = false;
            map<int,int> P;
            for(auto p : cdup->atom_parents) P[p]++;
            for(auto p : ptc) if (P[p]==2) {
                inside = true;
                if (P[p-1]<2 || P[p+1]<2) border = true;
            }
            if (border) return nullptr;
            if (inside) return cdup;
        }
        For(i, ptc.size()) ptc[i] = cdup->atom_parents[ptc[i]];
        if (cdup->type == "del") {
            delp.clear();
            for(int p : cdup->atom_parents) delp.insert(p);
            For(i, cdup->parent->atoms.size()) if (!delp.count(i)) ptc.push_back(i);
        }
        cdup = cdup->parent;

    }
    return nullptr;
}

double History::cherryness(const HAtom& a, const HAtom& b, int mode) {
    if (!mode) mode = cherry_mode;
    if (mode == KNOW_HOW) {
        HEvent* event = original->nth_from_end(events.size()-1);
        if (event->type == "del") event = resolve_deletion(event);
        if (event == nullptr) return 0.0;
        int x = -1, y = -2;
        For(i, event->atoms.size()) {
            if (event->atoms[i]==a) x = i;
            if (event->atoms[i]==b) y = i;
        }
        return 4.*((event->atom_parents[x] == event->atom_parents[y])
          &&((x > y) ^ event->compute_is_left()));
    }
    if (mode == NO_STRATEGY || cherry_forest == nullptr) return 1.0;
    if (mode == CHERRY_LEN) {
        double len = type_dna_length[a.atype()];
        return pow(cherry_forest->cherryness(a,b),
                   len/(len+1000));
    }
    return cherry_forest->cherryness(a,b);
}

void History::merge(const HAtom& a, const HAtom& b) {
    cherry_forest->merge(a,b);
}

void History::set_strategy(int strategy, Machine* machine) {
    switch(strategy) {
        case NO_STRATEGY:
        case CHERRY_NO:
        case CHERRY_TREE:
        case CHERRY_LEN:
        case KNOW_HOW:
            this->cherry_mode = strategy;
            assert(machine == nullptr);
            break;
        case SCORE_CL:
            this->cherry_mode = CHERRY_LEN;
            break;
        case SCORE_BAC_NC:
            this->cherry_mode = NO_STRATEGY;
            break;
        default:
            this->cherry_mode = DEFAULT_CHERRY;
    }
    this->strategy = strategy;
    this->machine = machine;
}

int History::get_history_score() {
    int score = 0;
    for (auto e : events) {
        auto type = e.second->type;
        if (type == "dup" || type == "dupi") {
            score += 1;
        } else if (type == "del") {
            score += 2;
        }
    }
    return score;
}
