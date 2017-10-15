#include"score.h"

typedef double (*scoring_function) (ScoringData* sd); 

double sc_seq_len(ScoringData* sd) {
    return log(SIZE(sd->last_event->atoms));
}
double sc_seq_num_types(ScoringData* sd) {
    return log(sd->num_types);
}
double sc_seq_dna_len(ScoringData* sd) {
    int len = 0;
    for(auto atom : sd->last_event->atoms)
        len += sd->history->type_dna_length[atom.atype()];
    return log(len);
}
double sc_post_bpc(ScoringData* sd) {
    return log(sd->post_bpc);
}
double sc_ev_len(ScoringData* sd) {
    return log(SIZE(sd->c.directions));
}
double sc_ev_dist(ScoringData* sd) {
    int e = min(sd->c.e1, max(sd->c.e2,sd->c.b2-1));
    int b = max(sd->c.b1+1, min(sd->c.e2, sd->c.b2+1));
    return log(b-e);
}
double sc_del_num(ScoringData* sd) {
    return SIZE(sd->deletions);
}
double sc_del_len(ScoringData* sd) {
    return log(1+SIZE(sd->c.directions) - SIZE(sd->atom_friends));
}
double sc_ev_sides(ScoringData* sd) {
    int dir = (sd->c.is_inv())?-1:1;
    int b1 = sd->c.b1, e1 = sd->c.e1+1;
    int b2 = sd->c.b2, e2 = sd->c.e2+dir;
    int size = SIZE(sd->last_event->atoms);
    int res = 0;
    if (b1 >= 0 && b2 >= 0 && b1<size && b2<size &&  
        sd->last_event->atoms[b1].type == dir*(sd->last_event->atoms[b2].type))
        res++;
    if (e1 >= 0 && e2 >= 0 && e1<size && e2<size &&  
        sd->last_event->atoms[e1].type == dir*(sd->last_event->atoms[e2].type))
        res++;
    return res;     
}
double sc_ev_prev_sp(ScoringData* sd) {
    return SIZE(sd->prev_sp)*0.5;
}
double sc_ev_post_sp(ScoringData* sd) {
    return SIZE(sd->post_sp)*0.5;
}
double sc_ev_prev_bp(ScoringData* sd) {
    return SIZE(sd->prev_bp)*0.5;
}
double sc_ev_post_bp(ScoringData* sd) {
    return SIZE(sd->post_bp)*0.5;
}
double sc_ev_dsig(ScoringData* sd) {
    if (!sd->c.is_inv()) return 0.0;
    set<int> prev_types, post_types;
    for(auto atom : sd->first_event->atoms) prev_types.insert(atom.type);
    for(auto atom : sd->last_event->atoms) post_types.insert(atom.type);
    int before = SIZE(prev_types) - sd->num_types;
    int diff = SIZE(post_types) - sd->num_types - before;
    return log(1.+diff/(1.+before));
}

double sc_avg_cherry(ScoringData* sd) {
    double sum = 0.;
    for(auto ap : sd->atom_friends) {
        sum += sd->history->cherryness(ap.first, ap.second, CHERRY_TREE);
    }
    return sum / double(SIZE(sd->atom_friends));
}
double sc_prod_cherry(ScoringData* sd) {
    double prod = 1.;
    for(auto ap : sd->atom_friends) {
        prod *= min(1.,sd->history->cherryness(ap.first, ap.second, CHERRY_TREE));
    }
    return prod;
}
double sc_len_cherry(ScoringData* sd) {
    double prod = 1.;
    for(auto ap : sd->atom_friends) {
        prod *= min(1.,sd->history->cherryness(ap.first, ap.second, CHERRY_LEN));
    }
    return prod;
}

vector<double> all_scores(History* h, const Candidate& c, HEvent* e) {
    vector<scoring_function> functions = {
        sc_seq_len,
        sc_seq_num_types,
        sc_post_bpc,

        sc_ev_len,
        sc_ev_dist,
        sc_del_num,
        sc_del_len,
        sc_ev_sides,
        sc_ev_prev_sp,
        sc_ev_post_sp,
        sc_ev_prev_bp,
        sc_ev_post_bp,
        sc_ev_dsig,

        sc_avg_cherry,
        sc_prod_cherry,
        sc_len_cherry,
    };
    vector<double> res(SIZE(functions));
    ScoringData* sd = new ScoringData(h,c,e);
    For(i, SIZE(functions)) {
        res[i] = functions[i](sd);
    }
    /*For(i, SIZE(functions)) {
        res.push_back(exp(res[i]));
        
        res.push_back(res[i]/(1+res[0]));
        res.push_back(res[i]/(1+res[1]));
        res.push_back(res[i]/(1+res[2]));
        res.push_back(res[i]/(1+res[3]));

        res.push_back(res[i]*res[0]);
        res.push_back(res[i]*res[1]);
        res.push_back(res[i]*res[2]);
        res.push_back(res[i]*res[3]);
    }*/

    delete sd;
    return res;
}

void compute_pairs(const vector<HAtom>& atoms, set<pii>& bp, set<pii>& sp) {
    int max_atom = 0, side = 0, dir = 0;
    for(auto atom : atoms) max_atom = max(max_atom, atom.atype());
    vi side_cnt[2] = {vi(max_atom + 1, 0), vi(max_atom + 1, 0)};
    vi side_last[2] = {vi(max_atom + 1, 0), vi(max_atom + 1, 0)};
    For(i, SIZE(atoms)+1) {
        int x = (i==0)?0:atoms[i-1].type;
        int y = (i==SIZE(atoms))?0:atoms[i].type;
        bp.insert({x, y});
        bp.insert({-y, -x});
        side = (x>=0);
        dir = (x<0)?-1:1;
        if (side_cnt[side][x*dir] == 0 || side_last[side][x*dir] != y*dir) {
            side_cnt[side][x*dir]++;
            side_last[side][x*dir] = y*dir;
        }
        side = (y<0);
        dir = (y<0)?-1:1;
        if (side_cnt[side][y*dir] == 0 || side_last[side][y*dir] != x*dir) {
            side_cnt[side][y*dir]++;
            side_last[side][y*dir] = x*dir;
        }
    }
    for(pii pp : bp) {
        if (side_cnt[pp.first>=0][abs(pp.first)] == 1 && side_cnt[pp.second<0][abs(pp.second)] == 1)
            sp.insert(pp);
    }
    /*cout << atoms << endl;
    cout << side_cnt[0] << endl << side_cnt[1] << endl;
    cout << "bp ";
    for(pii pp : bp) cout << pp.first << ":" << pp.second << "  ";
    cout << endl << "sp ";
    for(pii pp : sp) cout << pp.first << ":" << pp.second << "  ";
    cout << endl;*/
}

ScoringData::ScoringData(History* h, const Candidate& c, HEvent* e) : c(0,0,0,0) {
    this->history = h;
    this->c = c;
    h->rec_setup_scoring_data(c,e,this); 
    set<int> types;
    for(auto atom : last_event->atoms)
        types.insert(atom.atype());
    num_types = SIZE(types);
   
    compute_pairs(first_event->atoms, prev_bp, prev_sp);
    compute_pairs(last_event->atoms, post_bp, post_sp);
    prev_bpc = SIZE(prev_bp);
    post_bpc = SIZE(post_bp);
    set<pii> intersection_bp;
    for(pii bp : prev_bp) if (post_bp.count(bp)) intersection_bp.insert(bp);
    for(pii bp : intersection_bp) {
        prev_bp.erase(bp);
        post_bp.erase(bp);
    }
    set<pii> intersection_sp;
    for(pii sp : prev_sp) if (post_sp.count(sp)) intersection_sp.insert(sp);
    for(pii sp : intersection_sp) {
        prev_sp.erase(sp);
        post_sp.erase(sp);
    }

    int del_length = 0;
    for(pii d : deletions) del_length += d.second-d.first;
    assert(SIZE(c.directions) == SIZE(atom_friends) + del_length);
}

ScoringData::~ScoringData() {
    delete first_event;
    //delete second_event;
}
