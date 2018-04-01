#include"event.h"

Number Event::likelihood(){
    if (type=="root") 
        return Number(1.0);
    double lambda = history->parameters["likelihood_event_rate"]; 
    if (type=="sp" || type=="leaf"){
        return logNumber(-lambda * edge_time); 
    }
    double event_prob = lambda * exp(-lambda * edge_time);
    double begin_prob = 1.0/sequence_length;
    double p = 1.0/history->parameters["likelihood_mean_len"];
    double length_prob = p * pow(1.0-p, ev_len-1);
    p = 1.0/history->parameters["likelihood_mean_dist"];
    double dist_prob = p * pow(1.0-p, ev_dist);
    if (type=="dup"){
        double type_prob = (1.0 - history->parameters["likelihood_inversion_prob"])
                          *(1.0 - history->parameters["likelihood_deletion_prob"]);
        return Number(event_prob * type_prob * begin_prob * length_prob * dist_prob);
    }
    if (type=="dupi"){
        double type_prob = (history->parameters["likelihood_inversion_prob"])
                          *(1.0 - history->parameters["likelihood_deletion_prob"]);
        return Number(event_prob * type_prob * begin_prob * length_prob * dist_prob);
    }
    if (type=="del"){
        double type_prob = (1.0 - history->parameters["likelihood_deletion_prob"]);
        return Number(event_prob * type_prob * begin_prob * length_prob);
    }
    return Number(0.0);
}

void Event::read(istringstream& is){
    string parentname;
    is >> species >> name >> parentname >> event_time >> type;

    if (history->events.count(parentname)){
        parent = history->events[parentname];
    }else if (parentname=="root"){
        parent = NULL;
    }else {
        cerr << "Wrong parent name " << parent << endl;
    }
    history->events[name] = this;

    string token;
    while(is >> token){
        if (token[0] == '#') break;
        istringstream ss (token);
        int num;
        ss >> num;
        genes.push_back(num);
    }
    For(i, SIZE(genes)){
        int num;
        is >> num;
        gene_parents.push_back(num);
    }
}

void Event::write(ostream& os){
    os << species << " " << name << " " << ((parent==NULL)?"root":parent->name) 
       << " " << event_time << " " << type;
    for (auto g : genes) os << " " << g;
    os << " #";
    for (auto g : gene_parents) os << " " << g;
    os << endl;
}

void Event::writeshort(ostream& os){
    os << species << " " << name << " " << ((parent==NULL)?"root":parent->name) 
       << " " << event_time << " " << type;
    if (gene_parents[0] == -1)
        for (auto g : genes) os << " " << g;
    else{
        set<int> was;
        if (type=="del") {
            for (auto g : gene_parents) was.insert(g);
            For(i, SIZE(parent->genes)) if (!was.count(i))
                os << " " << parent->genes[i];
        } else for (auto g : gene_parents) {
            if (was.count(g)) os << " " << parent->genes[g];
            else was.insert(g);
        }
    }
    os << endl;
}

Event::Event(History* _history){
    this->history = _history;
    this->parent = NULL;
    int limit = SIZE(history->events)*2+5;
    do{ 
        this->name = "re" + to_string(rand()%limit);
        limit++;
    } while (history->events.count(this->name));
    history->events[this->name] = this;
}

