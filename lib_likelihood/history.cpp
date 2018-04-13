#include"history.h"
#include <iomanip>

namespace likelihood {
HistoryLikelihood::HistoryLikelihood(){
    kn_parameters = kn_atom_lists = kn_events = kn_trees = false;
    kn_alignments = 0;
}

void HistoryLikelihood::read_parameters(istream& is){
    string parameter_name;
    double parameter_value;
    while(is >> parameter_name >> parameter_value) {
        parameters[parameter_name] = parameter_value;
    }
    kn_parameters = true;
}

void HistoryLikelihood::set_parameters(double event_rate, double mean_dist, double mean_len,
    double inversion_prob, double deletion_prob) {
    parameters["likelihood_event_rate"] = event_rate;
    parameters["likelihood_mean_dist"] = mean_dist;
    parameters["likelihood_mean_len"] = mean_len;
    parameters["likelihood_inversion_prob"] = inversion_prob;
    parameters["likelihood_deletion_prob"] = deletion_prob;
    kn_parameters = true;
}

void HistoryLikelihood::read_atom_lists(istream& is){
    string species, id;
    int type, orientation, from, to;
    while(is >> species >> id >> type >> orientation >> from >> to) {
        if(atom_list.find(species) == atom_list.end())
            atom_list[species] = vector<atom_element>();
        atom_list[species].push_back(make_pair(id, type*orientation));
        existing_atoms.insert(type);
    }
    kn_atom_lists = true;
}

void HistoryLikelihood::read_atom_alignments(istream& is, int group_id){
    string name, word, line;

    while (is >> line) {
        if (line[0] == '>') {
            name = line.substr(1);

            atom_string[name] = "";
            atom_id[name] = group_id;
            group_atoms[group_id].push_back(name);
        } else {
            atom_string[name] += line;
        }
    }

    kn_alignments++;
}

void HistoryLikelihood::read_events(stringstream& is){
    string line;
    while(getline(is, line)){
        istringstream iss(line);
        new Event(this, iss);
    }
    kn_events = true;
}

inline bool compare_event_time(const Event* e1, const Event* e2){
    return e1->event_time < e2->event_time;
}

void HistoryLikelihood::write_events(ostream& os){
    os << std::setprecision(10);
    vector<Event*> ev;
    for(auto ep : events) ev.push_back(ep.second);
    sort(ev.begin(), ev.end(), compare_event_time);
    for(auto e : ev) e->writeshort(os);
}

void HistoryLikelihood::compute_events(){
    assert(kn_atom_lists && kn_alignments == SIZE(existing_atoms));
    for (auto x : atom_list) cerr << x.first << endl;

    if (SIZE(atom_list) != 1){
        cerr << "Neviem pocitat historiu viacerych druhov" << endl;
        return;
    }

    Event* state = new Event(this);
    state->load_leaf(*atom_list.begin());

    while(!state->final()){
        Event* e = state->find_prev_event();
        state = e;
    }
    state->gene_parents = vi(SIZE(state->genes), -1);
    double total_time = 0.0, wanted_time = branch_time;

    for(auto ep : events) {
        if (ep.second->type != "leaf")                      // leaf ma rovnake nody ako ten nad nim
            for (auto node : ep.second->atom_nodes) delete node;
        ep.second->atom_nodes.clear();

        if (ep.second->type == "root") total_time = -ep.second->event_time+0.01;
        ep.second->edge_time = ep.second->event_time;
    }
    for(auto ep : events) {
        if (ep.second->parent != NULL) ep.second->event_time = ep.second->parent->edge_time;
        else ep.second->event_time -= 0.01;
        ep.second->event_time = (ep.second->event_time*wanted_time/total_time) + wanted_time;
        if (ep.second->event_time < 1e-10) ep.second->event_time = 0;
        //cerr << "time " << ep.second->event_time << endl;
    }

    kn_events = true;
}

void HistoryLikelihood::clear_events(){
    for(auto ep : events) delete ep.second;
    events.clear();
    kn_events = false;
}


void HistoryLikelihood::build_trees(){
    assert(kn_atom_lists && kn_events);
    assert(kn_alignments == SIZE(existing_atoms));

    //zmaz pripadne stare stromceky
    for(auto ap : group_trees) delete ap.second;
    group_trees.clear();

    //vyrob atomove stromceky
    for(auto ep : events) {
        //zmaz pripadne stare stromy
            ep.second->atom_nodes.clear();
        For(i,ep.second->genes.size()) {
            ep.second->atom_nodes.push_back(new Node());
            ep.second->atom_nodes[ep.second->atom_nodes.size()-1]->event = ep.second;
        }
        if(ep.second->parent == NULL) ep.second->edge_time = 0;
        else ep.second->edge_time = ep.second->event_time - ep.second->parent->event_time;
    }
    for(auto ep : events) {
        Event* e = ep.second;
        For(i,e->atom_nodes.size()) {
            if(e->gene_parents[i]==-1) continue;
            e->atom_nodes[i]->edge_time = e->edge_time;
            e->atom_nodes[i]->parent = e->parent->atom_nodes[e->gene_parents[i]];
            if(e->parent->atom_nodes[e->gene_parents[i]]->left == NULL)
                e->parent->atom_nodes[e->gene_parents[i]]->left = e->atom_nodes[i];
            else
                e->parent->atom_nodes[e->gene_parents[i]]->right = e->atom_nodes[i];
        }
        if (e->type == "leaf") {
            For(i, e->atom_nodes.size()){
                e->atom_nodes[i]->data = atom_string[atom_list[e->species][i].first];
            }
        }
    }
    for(auto ep : events) {
        Event* e = ep.second;
        if(e->type != "root") continue;
        For(i,e->atom_nodes.size()) {
            int typ = abs(e->genes[i]);
            AtomTree* at = new AtomTree();
            at->root = e->atom_nodes[i];
            if (existing_atoms.count(typ))
                at->atom_length = atom_string[group_atoms[typ][0]].length();
            else
                at->atom_length = 0;
            group_trees[typ] = at;
        }
    }
    //porataj ev_len, ev_dist a sequence_length
    map<int, int> parent_last;
    for (auto ep : events){
        Event* e = ep.second;
        if (e->type == "root" || e->type == "sp" || e->type == "leaf") continue;
        int parentsum = 0;
        for(auto x : e->parent->genes)
            parentsum+=group_trees[abs(x)]->atom_length;
        e->sequence_length = parentsum;

        parent_last.clear();
        int sum = 0, minsum = 2*parentsum, maxsum = 0;
        for (auto x : e->gene_parents){
            if (parent_last.count(x)) {
                maxsum = max(maxsum, parent_last[x]);
                minsum = min(minsum, sum);
            }
            sum += group_trees[abs(e->parent->genes[x])]->atom_length;
            parent_last[x] = sum;
        }
        if (e->type == "del"){
            e->ev_len = parentsum-sum;
            e->ev_dist = 0;
        }else{
            e->ev_len = sum-parentsum;
            e->ev_dist = minsum - maxsum;
        }
        //cerr << e->name << e->type << " " << e->sequence_length << " " << e->ev_len << " " << e->ev_dist << endl;
    }

    kn_trees = true;
}

Number HistoryLikelihood::likelihood_all(){
    assert(kn_trees && kn_parameters);

    Number total = Number(1.0);
    // zrataj likelihood kazdeho atomoveho stromu
    if (felsenstein_on){
        for (auto t : group_trees){
            Number x = t.second->likelihood();
            //x.print(stderr);
            total *= x;
        }
        //fprintf(stderr, ".. zatial: "); total.print(stderr);
    }
    // zrataj likelihood eventov
    for (auto e : events){
        Number x = e.second->likelihood();
        //fprintf(stderr, "%s: ", e.second->name.c_str());
        //x.print(stderr);
        total *= x;
    }
    return total;
}


void HistoryLikelihood::add_to_old(Event* e, double koef) {
    if (e->type=="dup" || e->type=="dupi"){
        old_seqences[e->parent->genes]+=koef;
        old_events[e->get_change()]+=koef;
    }
}
}
