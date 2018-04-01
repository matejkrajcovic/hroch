#include"event.h"
#include"atomtree.h"

void Event::load_leaf(const pair<string, vector<atom_element> >& atom_list){
    this->parent = NULL;
    this->species = atom_list.first;
    this->type = "leaf";
    this->event_time = 0.0;

    this->genes.clear();
    for(auto ae : atom_list.second)
        this->genes.push_back(ae.second);

    this->gene_parents.resize(SIZE(genes));
    For(i, SIZE(genes)) this->gene_parents[i] = i;
    
    if (history->scoring_similar) {
        For(g, SIZE(genes)) {
            Node *node = new Node();
            node->data = history->atom_string[history->atom_list[this->species][g].first];
            //cerr << "leaf " << g << " " << genes[g] << " " << SIZE(node->data) << endl;
            For(c, 4){
                node->allprobs[c].resize(SIZE(node->data));
                For(i, SIZE(node->allprobs[c])) 
                    node->allprobs[c][i] = (node->data[i]=='N' || node->data[i]=='-' || 
                        node->data[i]==bases[c])?1.:0.;
            }
            node->data.clear();
            this->atom_nodes.push_back(node);
        }
    }
}

bool Event::final(){
    set<int> types;
    for(auto x : genes) types.insert(abs(x));
    return (SIZE(types) == SIZE(genes));
}

// pozor, zalezi na poradi jednotlivych volani tejto funckie
void Event::add_edge(int fromk, int fromi, int fromj, int tok, int toi, int toj, double val){
    val *= S[fromk][fromi][fromj];
    if (val < 1e-9) return;
    P[tok][toi][toj].push_back(SIZE(E));
    E.push_back(make_pair(trint(fromk, fromi, fromj), val));
    S[tok][toi][toj] += val; 
}

trint Event::path_previous(trint x){
    double r = random_double(0.0, S[x.a][x.b][x.c]);
    double sum = 0;
    for (auto v : P[x.a][x.b][x.c]){
        sum += E[v].second;
        if (r<sum) return E[v].first;
    }
    assert(0);
}


inline void genes_cpy(Event* from, Event* to, int b, int e){
    to->genes.insert(to->genes.end(), from->genes.begin()+b, from->genes.begin()+e);
}
    
void Event::compute_prevgenes(const vtri& path, Event* e, vector<piii>& deletions, vector<int>& data){ //{{{
    int pfb = path.front().b, pfc = path.front().c,
        pbb = path.back().b+1, pbc = path.back().c+1;
    int siga = (path.front().a==0)?1:-1;
    int posa = (siga>0)?0:3, posb = pfb-1, posc = pfc-1*siga;
    int adi1, adi2, adi3, adi4, adi5; // dlzky usekov v e; adi2*adi4 = 0
    if (siga < 0){
        pbc--; pfc++;
        swap(pfc, pbc);
    }

    genes_cpy(this, e, 0, min(pfb, pfc));
    adi1 = SIZE(e->genes);
    if(pfb >= pfc) genes_cpy(this, e, pbc, pfb); // ked duplikujem zprava dolava
    adi2 = SIZE(e->genes)-adi1;

    deletions.clear();
    for(auto x : path){
        if (x.a % 3) {
            if (x.a != posa)
                deletions.push_back(piii(pii(((x.a%3==1)?x.c+(siga>0):x.b+1), 
                                         SIZE(e->genes)), 0));
            deletions.back().second++;
        }
        if (x.b != posb) e->genes.push_back(genes[x.b]);
        else if (x.c != posc) e->genes.push_back(genes[x.c]*siga);
        posa = x.a; posb = x.b; posc = x.c;
    }
    adi3 = SIZE(e->genes)-adi2-adi1;
    
    if(pfc >= pfb) genes_cpy(this, e, pbb, pfc); // ked duplikujem zlava doprava
    adi4 = SIZE(e->genes)-adi3-adi2-adi1;
    genes_cpy(this, e, max(pbb, pbc), SIZE(genes));
    adi5 = SIZE(e->genes)-adi4-adi3-adi2-adi1;
    data = vector<int> {siga, adi1, adi2, adi3, adi4, adi5, pfb, pbb, pfc, pbc};
} //}}}

Event* Event::event_from_path(const vtri& path){ //{{{
    Event* e = new Event(history);
    e->parent = NULL;
    e->species = this->species;
    e->type = "root";    
    
    // e je udalost pred duplikaciou, teraz jej vyrobime genes
    
    vector<piii> deletions; // suradnice su podla ((this.zac, e.zac) len)
    vector<int> data; // siga, adi1..adi5, pfb, pbb, pfc, pbc
    compute_prevgenes(path, e, deletions, data);
    int &siga=data[0], &pfb = data[6], &pbb = data[7], &pfc = data[8];
    int &adi1 = data[1], &adi2 = data[2], &adi3 = data[3], &adi4 = data[4], &adi5 = data[5];

    double gap_time = 0.01;
    if (history->scoring_similar){
        vector<pair<Node*, Node*> > nodepairs;
        for(auto x : path) if (x.a%3==0){ 
            assert(SIZE(atom_nodes[x.b]->allprobs[0])==SIZE(atom_nodes[x.c]->allprobs[0]));
            nodepairs.push_back(make_pair(atom_nodes[x.b], atom_nodes[x.c])); 
        }
        gap_time = gap_time*0.2 + best_time(nodepairs)*0.8;
    }

    // vsunieme medzi jednotlive duplikacie
    Event* laste = this;
    sort(deletions.begin(), deletions.end(), greater<piii>());
    for(auto x : deletions){
        Event* nevent = new Event(history);
        nevent->species = laste->species;
        nevent->parent = NULL;
        nevent->event_time = laste->event_time-gap_time/(SIZE(deletions)+1.);
        genes_cpy(laste, nevent, 0, x.first.first); 
        if (siga>0 || (x.first.first >= pfb && x.first.first <= pbb))
            genes_cpy(e, nevent, x.first.second, x.first.second + x.second); 
        else
            For(i, x.second) nevent->genes.push_back(-e->genes[x.first.second+x.second-i-1]);
        genes_cpy(laste, nevent, x.first.first, SIZE(laste->genes)); 

        laste->gene_parents.resize(SIZE(laste->genes));
        For(i, SIZE(laste->genes)) laste->gene_parents[i] = i + (i>=x.first.first)*x.second;

        laste->type = "del";
        laste->parent = nevent;
        laste = nevent;
    }
   
    e->event_time = laste->event_time-gap_time/(SIZE(deletions)+1.);
    laste->parent = e;
    if (siga>0) laste->type = "dup";
    else laste->type = "dupi";

    // popocitaj parentov
    laste->gene_parents.resize(SIZE(laste->genes));
    For(i, adi1) laste->gene_parents[i] = i;
    if (siga>0 || pfb < pfc) For(i, adi3) laste->gene_parents[adi1+i] = adi1+adi2+i;
    else For(i, adi3) laste->gene_parents[adi1+adi3-i-1] = adi1+adi2+i;
    For(i, adi2) laste->gene_parents[adi1+adi3+i] = adi1+i;
    For(i, adi4) laste->gene_parents[adi1+adi3+adi2+i] = adi1+adi2+adi3+i;
    if (siga>0 || pfc < pfb) For(i, adi3) laste->gene_parents[adi1+adi3+adi2+adi4+i] = adi1+adi2+i; 
    else For(i, adi3) laste->gene_parents[adi1+adi3+adi2+adi4+adi3-i-1] = adi1+adi2+i; 
    For(i, adi5) laste->gene_parents[adi1+adi3+adi2+adi4+adi3+i] = adi1+adi2+adi3+adi4+i;

    // popocitaj allprobs, ked je to zapnute
    if (history->scoring_similar){
        For(i, SIZE(e->genes)) e->atom_nodes.push_back(new Node());
        For(i, SIZE(genes)) {
            int pos = i;
            Event* event = this;
            while(event!=e) {
                pos = event->gene_parents[pos];
                event = event->parent;
            }
            atom_nodes[i]->parent = e->atom_nodes[pos];
            if (atom_nodes[i]->parent->left == NULL)
                atom_nodes[i]->parent->left = atom_nodes[i];
            else
                atom_nodes[i]->parent->right = atom_nodes[i];
        }
        For(g, SIZE(e->genes)){
            Node* node = e->atom_nodes[g];
            assert(node->left!=NULL);
            For(c, 4){
                node->allprobs[c].clear();
                node->allprobs[c].resize(SIZE(node->left->allprobs[c]));
                For(i, SIZE(node->allprobs[c])) {
                    if (node->right==NULL) {
                        For(j, 4)
                            node->allprobs[c][i] += change_probability(c, j, gap_time) *
                                node->left->allprobs[j][i];
                    } else {
                        For(j, 4)
                            node->allprobs[c][i] += change_probability(c, j, gap_time) *
                                (node->left->allprobs[j][i]+node->right->allprobs[j][i]) * 0.5;
                    }
                }
            }
            node->data.clear();
        }
    }
    return e;
} //}}}
    
void Event::compute_graph(){ //{{{
    // Mame graf, ktoreho vrcholy su usporiadane v 6 stvorcoch n x n, ku ktorym je 
    // pridana siedma vrstva s vrcholmi beg a end.
    // Jednotlive stovrce maju postupne funkciu S, D, R, -S, -D, -R
    
    double c_dup = 1.0,
           c_dupi = 1.0,
           c_dupm = 1.8,
           c_del = 0.22,
           c_delm = 0.8;

    int n = SIZE(genes);
  
    S = vector<vvdo>(6, vvdo(n, vdo(n, 0.0)));
    S.push_back(vvdo(1,vdo(2, 0.0)));
    P = vector<vvvi>(6, vvvi(n, vvi(n, vi(0))));
    P.push_back(vvvi(1, vvi(2, vi(0))));
    E.clear();
    
    S[beg.a][beg.b][beg.c] = 1;
    For(i, n) For(j, n) {                                                       // duplikacia s deleciami
        if (i!=j && genes[i] == genes[j]){
            add_edge(beg.a,beg.b,beg.c, 0,i,j, c_dup);
            if (i>0 && j>0) {
                add_edge(1,i-1,j-1, 0,i,j, c_del);
                add_edge(2,i-1,j-1, 0,i,j, c_del);
                if(genes[i-1] == genes[j-1])
                    add_edge(0,i-1,j-1, 0,i,j, c_dupm);
            }
            if (i<n-1) add_edge(0,i,j, 1,i+1,j, c_del);
            if (j<n-1) add_edge(0,i,j, 2,i,j+1, c_del);
            add_edge(0,i,j, end.a,end.b,end.c, c_dup);
        }
        if (i>0 && i!=j) {
            add_edge(2,i-1,j, 1,i,j, c_del*c_del);
            add_edge(1,i-1,j, 1,i,j, c_delm);
        }
        if (j>0 && i!=j) {
            add_edge(1,i,j-1, 2,i,j, c_del*c_del);
            add_edge(2,i,j-1, 2,i,j, c_delm);
        }
    }
    For(i, n) for(int j = n-1; j>=0; --j) {                                     // inverzna duplikacia s deleciami
        if (i!=j && genes[i] == -genes[j]){
            add_edge(beg.a,beg.b,beg.c, 3,i,j, c_dupi);
            if (i>0 && j<n-1 && i!=j+1) {
                add_edge(4,i-1,j+1, 3,i,j, c_del);
                add_edge(5,i-1,j+1, 3,i,j, c_del);
                if(genes[i-1] == -genes[j+1])
                    add_edge(3,i-1,j+1, 3,i,j, c_dupm);
            }
            if (i<n-1 && i+1!=j) add_edge(3,i,j, 4,i+1,j, c_del);
            if (j>0   && i+1!=j) add_edge(3,i,j, 5,i,j-1, c_del);
            add_edge(3,i,j, end.a,end.b,end.c, c_dupi);
        }
        if (i>0 && i!=j) {
            add_edge(5,i-1,j, 4,i,j, c_del*c_del);
            add_edge(4,i-1,j, 4,i,j, c_delm);
        }
        if (j<n-1 && i!=j) {
            add_edge(4,i,j+1, 5,i,j, c_del*c_del);
            add_edge(5,i,j+1, 5,i,j, c_delm);
        }
    }

    /*For(i, n) {
        For(j, n)
            if (S[0][i][j]>1e-5) printf("%5.2lf ", S[0][i][j]);
            else printf("  .   ");
        printf("   |   ");
        For(j, n)
            if (S[1][i][j]>1e-5) printf("%5.2lf ", S[1][i][j]);
            else printf("  .   ");
        printf("   |   ");
        For(j, n)
            if (S[2][i][j]>1e-5) printf("%5.2lf ", S[2][i][j]);
            else printf("  .   ");
    printf("\n"); }
    printf("\n");
    For(i, n) {
        For(j, n)
            if (S[3][i][j]>1e-5) printf("%5.2lf ", S[3][i][j]);
            else printf("  .   ");
        printf("   |   ");
        For(j, n)
            if (S[4][i][j]>1e-5) printf("%5.2lf ", S[4][i][j]);
            else printf("  .   ");
        printf("   |   ");
        For(j, n)
            if (S[5][i][j]>1e-5) printf("%5.2lf ", S[5][i][j]);
            else printf("  .   ");
    printf("\n"); }*/
} //}}}

vtri Event::apply_event(const vi& event){
    //int k = 0;
    //for(auto x : genes) cout << x << ((++k%10)?" ":" | "); cout << endl;

    int pos1 = event[SIZE(event)-2];
    int pos2 = event[SIZE(event)-1];
    int dir = 1, offs = 0;
    if (pos2<0) {pos2 = -pos2; dir = -1; offs = 3;}
    vtri path;
    For(i, SIZE(event)-2){
        if (pos1 >= SIZE(genes) || pos2 >= SIZE(genes) || pos2 < 0){
            path.clear();
            return path;
        }else if (genes[pos1] == dir*genes[pos2]){
            path.push_back(trint(offs, pos1, pos2));
            pos1++; pos2+=dir;
        }else if (genes[pos1] == event[i] && i!=0 && i!=SIZE(event)-3){
            path.push_back(trint(offs+1, pos1, pos2-dir));
            pos1++;
        }else if (dir*genes[pos2] == event[i] && i!=0 && i!=SIZE(event)-3){
            path.push_back(trint(offs+2, pos1-1, pos2));
            pos2+=dir;
        }else{
            path.clear();
            return path;
        }
    }
    return path;
}

Event* Event::find_prev_event(){ //{{{
    if (type=="leaf") {
        Event* e = new Event(history);
        e->species = this->species;
        e->event_time = this->event_time;
        e->genes = this->genes;
        e->atom_nodes = this->atom_nodes;
        this->parent = e;
        return e;
    }

    compute_graph();

    map<vtri, double> ways;
    double waysum = 0.0;

    int iterations = 0;
    while(iterations++ <= SIZE(genes)*5){
        vtri path;
        trint x = end;
        while(x != beg){
            path.push_back(x);
            x = path_previous(x);
        }
        reverse(path.begin(), path.end());
        path.pop_back();

        assert(SIZE(path)>0);
        
        int f1 = path.front().b, b1 = path.back().b,
            f2 = min(path.front().c, path.back().c),
            b2 = max(path.front().c, path.back().c);
        if (min(b1, b2) >= max(f1, f2)) continue;
        if (ways.count(path)) continue;
        //cerr << f1 << " " << b1 << " -> " << f2 << " " << b2 << endl;

        ways[path] = path_score(path);
        if (SIZE(ways)>=SIZE(genes)/2) break;
    }
    if (history->scoring_oldhistory){
        double best = 0.0;
        for(auto x : history->old_events) 
            best = max(best, x.second);

        int cnt = 0;
        for(auto x : history->old_events) {
            if (x.second < 0.5*best) continue;
            cnt++;
            vtri path = apply_event(x.first);
            if (SIZE(path)){
                //cerr << "pth "; for(auto y : x.first) cerr << y << " "; cerr << endl;
                if (!ways.count(path)) ways[path] = path_score(path);
            }
        }
        //cout << cnt << " " << SIZE(genes) << endl;
    }
    assert(SIZE(ways));
    for (auto x : ways){
        fprintf(stderr, "score %10.2lf            %3d %3d -> %3d %3d\n", x.second,
               x.first.front().b, x.first.back().b,
               x.first.front().c, x.first.back().c);
        waysum += x.second;
    }

    double x = random_double(0.0, waysum);
    cerr << "chosen " << x << endl;
    waysum = 0.0;
    for(auto pp : ways) {
        waysum+=pp.second;
        if (x < waysum) return event_from_path(pp.first);
    }
    assert(0);
    return NULL;
} //}}}

vi Event::get_change(){
    if (type=="leaf" || type=="root" || type=="sp" || type=="del")
        return vi(0);
    map<int, int> was;
    int fi=-1, se=0; 
    vi result;
    For(i, SIZE(gene_parents)){
        if (was.count(gene_parents[i])) {
            result.push_back(parent->genes[gene_parents[i]]);
            if (fi<0){
                fi = was[gene_parents[i]];
                se = i;
                if (gene_parents[fi]!=fi) swap(fi, se);
                if (genes[fi] == -genes[se]) se = -se;
            }
        }
        was[gene_parents[i]] = i;
    }
    result.push_back(fi);
    result.push_back(se);
    return result;
}
