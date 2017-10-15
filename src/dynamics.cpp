#include"dynamics.h"

#define BASE 0
#define BASEI 4
#define LEV 1
#define LEVY 2
#define LEVX 3

Dynamics::Dynamics(History* history, HEvent* event) {
    this->event = event;
    for(auto a : event->atoms) atoms.push_back(a.type);
    n = SIZE(atoms);
    cherryness = vector<vvdo>(n+2, vvdo(n+2, vdo(8, 0.)));
    For(i, n) For(j, n) {
        if (i!=j && atoms[i]==atoms[j]) cherryness[i+1][j+1][BASE+LEV] = 
            history->cherryness(event->atoms[i],event->atoms[j]);
        if (i!=j+1 && atoms[i]==-atoms[j]) cherryness[i+1][j+1][BASEI+LEV] = 
            history->cherryness(event->atoms[i],event->atoms[j]);
        For(k,2) 
            cherryness[i+1][j+1][BASE+k+2] = cherryness[i+1][j+1][BASEI+k+2] = 1.;
    }
}

//                 B   D DY DX BI  I IY IX
const int dx[] = { 0, -1,-1, 0, 0,-1,-1, 0};
const int dy[] = { 0, -1, 0,-1, 0, 1, 0, 1};

// Use extreme caution in next two functions, 
// they are weird to be as fast as possible.
void Dynamics::compute_graph(double dupm, double delm, double deli) {
    this->dupm = dupm;
    this->delm = delm;
    this->deli = deli;
    this->deli3 = deli*deli*deli;
    all_endpoints.clear();
    mass = vector<vvdo>(n+2, vvdo(n+2, vdo(8, 0.0)));
    const double multi[4][4] = {
        {0.,0.,0.,0.},  // base
        {startp,dupm,dupm,dupm},   
        {0.,deli*delm,delm,deli3*delm},
        {0.,deli*delm,deli3*delm,delm},
    };
    for(int i = 1; i<=n; ++i) {
        // duplications
        for(int j = 1; j<=n; ++j) if (i!=j) {
            mass[i+dx[BASE+LEV]][j+dy[BASE+LEV]][BASE] = 1.0;
            for(int k = 1; k<4; ++k)
                for(int l = 0; l<4; ++l)
                    mass[i][j][BASE+k] += mass[i+dx[BASE+k]][j+dy[BASE+k]][BASE+l] * 
                        multi[k][l] * cherryness[i][j][BASE+k];
            if(i!=j && atoms[i-1] == atoms[j-1]) 
                all_endpoints.push_back({mass[i][j][BASE+LEV],trint(i,j,BASE+LEV)});
        }
        // inverse duplications
        for(int j = n; j>0; --j) if (i!=j) {
            mass[i+dx[BASEI+LEV]][j+dy[BASEI+LEV]][BASEI] = 1.0;
            for(int k = 1; k<4; ++k)
                for(int l = 0; l<4; ++l)
                    mass[i][j][BASEI+k] += mass[i+dx[BASEI+k]][j+dy[BASEI+k]][BASEI+l] 
                        * multi[k][l] * cherryness[i][j][BASEI+k];
            if(atoms[i-1] == -atoms[j-1]) 
                all_endpoints.push_back({mass[i][j][BASEI+LEV],trint(i,j,BASEI+LEV)});
        }
    }
    // optimalisation
    sort(all_endpoints.begin(), all_endpoints.end());
    reverse(all_endpoints.begin(), all_endpoints.end());
    sum_endpoints = 0.;
    for(auto ep : all_endpoints) sum_endpoints += ep.first;
    /*csum_ep = vector<double>(SIZE(all_endpoints)+1,0.);
    For(i, SIZE(all_endpoints))
        csum_ep[i+1] = csum_ep[i]+all_endpoints[i].first;*/
}

Candidate Dynamics::get_candidate(bool maximal) {
    assert(SIZE(all_endpoints));
    const double multi[4][4] = {
        {0.,0.,0.,0.},
        {startp,dupm,dupm,dupm},
        {0.,deli*delm,delm,deli3*delm},
        {0.,deli*delm,deli3*delm,delm},
    };

    int pos = 0;
    if (!maximal && sum_endpoints > 0) {
       double r = random_double(0., sum_endpoints);
       // its faster this way
       for(pos = 0; r>all_endpoints[pos].first; pos++)
           r-=all_endpoints[pos].first;
       //pos = upper_bound(csum_ep.begin(), csum_ep.end(), r) - csum_ep.begin() - 1;
    } else {
        pos = 0;
    }
    
    trint xyz = all_endpoints[pos].second;
    int x = xyz.a, y = xyz.b, z = xyz.c;
    int end1 = x, end2 = y;
    vector<int> directions;
    
    while(z%4>0) {
        if (mass[x][y][z] < epsilon) return Candidate(-1,0,-1,0);
        double r = random_double(0., mass[x][y][z]);
        if (maximal) {
            double max_val = -1;
            double sum_val = epsilon, val = 0.0;
            For(k,4) {
                int z2 = (z/4)*4 + k;
                val = mass[x+dx[z]][y+dy[z]][z2] * multi[z%4][k] * cherryness[x][y][z];
                if (max_val < val) {
                    max_val = val;
                    r = sum_val;
                }
                sum_val += val;
            }
        }

        For(k,4) {
            int z2 = (z/4)*4 + k;
            if ((r -= mass[x+dx[z]][y+dy[z]][z2] * multi[z%4][k] * cherryness[x][y][z]) < 0.) {
                directions.push_back(z%4-1); 
                x += dx[z];
                y += dy[z];
                z = z2;
                break;
            }
        }
    }
    reverse(directions.begin(), directions.end());
    return Candidate(x-1,end1-1,y-1,end2-1,directions); 
}

void compute_A0A1P1(const Candidate& C, const vector<HAtom>& atoms, 
    vector<HAtom>& A0, vector<HAtom>& A1, vector<int>&P1, vector<pii>& deletions) { // {{{
    vector<HAtom> merged[3];
    deletions.clear();
    int dir[3] = {0,1,C.is_inv()?-1:1};
    int p[3] = {0, C.b1, C.b2}, lastd = 0, lastlen = 0;
    for(auto d : C.directions) {
        p[0]++;
        if (d) {
            p[d] += dir[d];
            merged[d].push_back(atoms[p[d]]);
            merged[0].push_back(atoms[p[d]]);
            merged[0].back().type*=dir[d];
            merged[3-d].push_back(HAtom(atoms[p[d]].type*dir[2]));
        } else {
            p[1] += dir[1];
            p[2] += dir[2];
            merged[1].push_back(atoms[p[1]]);
            merged[2].push_back(atoms[p[2]]);
            merged[0].push_back(atoms[p[1]]);
            merged[0].back().add_ids(atoms[p[2]].get_ids());
        }
        lastlen++;
        if (lastd != d) {
            if (lastd) {
                int a = ((p[lastd] < p[3-lastd])?-1:1)*(p[0]-lastlen-1);
                int b = lastlen*dir[3-lastd];
                deletions.push_back({a,b});
            }
            lastlen = 0;
            lastd = d;
        }
    }
    if (C.is_inv()) reverse(merged[2].begin(), merged[2].end());

    int oa = min(C.b1+1,min(C.b2+1,C.e2));
    int ob = min(C.e1+1,max(C.e2+1,C.b2));
    int oc = max(C.b1+1,min(C.b2+1,C.e2));
    int od = max(C.e1+1,max(C.e2+1,C.b2));
    int oe = SIZE(atoms);
    int ooam = oa+(oc-ob)*C.is1_right();
    int doom = oa+oc+SIZE(merged[0])-ob;
    int x1 = C.is1_left()?1:2;
    for(auto& d : deletions) {
        //cout << d.first << " " << d.second << endl;
        if (d.second < 0) {
            d.second *= -1;
            d.first = (SIZE(merged[0]) - abs(d.first) - d.second)*sign(d.first);
        }
        d.first = ((d.first<0)?doom-d.first:oa+d.first);
        d.second = d.first + d.second;
    }
    P1.reserve(oe);
    
    // first part
    A0.clear(), A1.clear();
    For(i, oa) P1.push_back(SIZE(A0)+i);
    A0.insert(A0.end(), atoms.begin(), atoms.begin()+oa);
    A1.insert(A1.end(), atoms.begin(), atoms.begin()+oa);

    // merge1
    if (C.is1_left()) A0.insert(A0.end(), merged[0].begin(), merged[0].end());
    if (C.is1_left() || !C.is_inv()) For(i, SIZE(merged[0])) P1.push_back(ooam+i);
    else For(i, SIZE(merged[0])) P1.push_back(ooam+SIZE(merged[0])-1-i);
    A1.insert(A1.end(), merged[x1].begin(), merged[x1].end());
    
    For(i, oc-ob) P1.push_back(SIZE(A0)+i);
    A0.insert(A0.end(), atoms.begin()+ob, atoms.begin()+oc);
    A1.insert(A1.end(), atoms.begin()+ob, atoms.begin()+oc);
    
    // merge2
    if (C.is1_right()) A0.insert(A0.end(), merged[0].begin(), merged[0].end());
    if (C.is1_right() || !C.is_inv()) For(i, SIZE(merged[0])) P1.push_back(ooam+i);
    else For(i, SIZE(merged[0])) P1.push_back(ooam+SIZE(merged[0])-1-i);
    A1.insert(A1.end(), merged[3-x1].begin(), merged[3-x1].end());
    
    // last part
    For(i, oe-od) P1.push_back(SIZE(A0)+i);
    A0.insert(A0.end(), atoms.begin()+od, atoms.begin()+oe);
    A1.insert(A1.end(), atoms.begin()+od, atoms.begin()+oe);
} // }}}

HEvent* compute_next_AP(const pii& deletion, const string& species, const string& name,
        HEvent* parent, vector<HAtom>& A, vector<int>& P) {
    HEvent* event = new HEvent(species, name, "");
    event->parent = parent;
    event->atoms = A;
    event->atom_parents = P;
    A.clear();
    P.clear();
    For(i, SIZE(event->atoms)) if (i<deletion.first || i>=deletion.second) {
        A.push_back(event->atoms[i]);
        P.push_back(i);
    }
    return event;
}

void History::rec_compute_parent(const Candidate& C, HEvent* event) {
    vector<HAtom> A0, A;
    vector<int> P;
    vector<pii> deletions;

    compute_A0A1P1(C, event->atoms, A0,A,P, deletions);
    sort(deletions.begin(), deletions.end(), greater<pii>());
    
    HEvent* current = new HEvent(event->species, gen_event_name(), "");
    current->atoms = A0;
    for(auto d : deletions) { 
        current = compute_next_AP(d, event->species, gen_event_name(), 
                                  current, A, P);
        current->type = (current->parent->type=="")?(C.is_inv()?"dupi":"dup"):"del";
    }

    assert(event->atoms == A);
    event->atom_parents = P;
    event->parent = current;
    event->type = (event->parent->type=="")?(C.is_inv()?"dupi":"dup"):"del";

    current = event;
    while(current->parent != nullptr) {
        current->parent->compute_atom_ids(current);
        current = current->parent;
    }
}

HEvent* History::rec_see_event(const Candidate& C, HEvent* event) {
    vector<HAtom> A0, A;
    vector<int> P;
    vector<pii> deletions;
    compute_A0A1P1(C, event->atoms, A0,A,P, deletions);
    HEvent* e = new HEvent(event->species, "test", "dup");
    if (C.is_inv()) e->type += "i";
    e->atoms = A;
    e->atom_parents = P;
    return e;
}

void History::rec_setup_scoring_data(const Candidate& C, HEvent* event, ScoringData* sd) {
    vector<HAtom> A0, A;
    vector<int> P;
    compute_A0A1P1(C, event->atoms, A0,A,P, sd->deletions);
    sd->first_event = new HEvent(event->species, "test", "");
    sd->first_event->atoms = A0;
    // Uncomment if you need second event too:
    /*sd->second_event = new HEvent(event->species, "test", "dup");
    if (C.is_inv()) sd->second_event->type += "i";
    sd->second_event->atoms = A;
    sd->second_event->atom_parents = P;*/
    sd->last_event = event; 
    
    int dir[3] = {0,1,C.is_inv()?-1:1};
    int p[3] = {0, C.b1, C.b2};
    for(auto d : C.directions) {
        if (d) {
            p[d] += dir[d];
        } else {
            p[1] += dir[1];
            p[2] += dir[2];
            sd->atom_friends.push_back({event->atoms[p[1]], event->atoms[p[2]]});
        }
    }    
}

void History::rec_merge_candidate(const Candidate& C, HEvent* event) {
    int dir[3] = {0,1,C.is_inv()?-1:1};
    int p[3] = {0, C.b1, C.b2};
    for(auto d : C.directions) {
        if (d) {
            p[d] += dir[d];
        } else {
            p[1] += dir[1];
            p[2] += dir[2];
            merge(event->atoms[p[1]], event->atoms[p[2]]);
        }
    }    
}
