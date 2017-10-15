#include"cherry.h"

bool CherryTree::has_cherry(const HAtom& a, const HAtom& b) {
    int x = id[a], y = id[b];
    return edges[x] == edges[y];
}

void CherryTree::merge(const HAtom& a, const HAtom& b) {
    int x = id[a], y = id[b];
    int z = edges[x][0];
    auto it = find(edges[z].begin(), edges[z].end(), x);
    edges[z].erase(it);
    it = find(edges[z].begin(), edges[z].end(), y);
    edges[z].erase(it);

    HAtom c = a;
    c.add_ids(b.get_ids());
    id[c] = z;
}

CherryTree::CherryTree(const CherryTree& ct) {
    this->id = ct.id;
    this->edges = ct.edges;
    this->atype = ct.atype;
    this->probability = ct.probability;
}

CherryTree::CherryTree(double probability, int atype, map<string, int>& ids, const string& word) {
    this->probability = probability;
    this->atype = atype;
    this->parse(ids, word, 0, SIZE(word));
}

int CherryTree::parse(map<string,int>& ids, const string& word, int from, int to) {
    int me = SIZE(edges);
    edges.push_back(vi());

    if (word[from] == '(') {
        int pos = from+1, depth = 0;
        while(pos < to) {
            if (word[pos] == '(') depth++;
            if (word[pos] == ')') depth--;
            if ((depth == 0 && word[pos] == ',') || depth < 0) {
                int son = parse(ids, word, from+1, pos); 
                edges[me].push_back(son);
                edges[son].push_back(me);
                from = pos;
            }
            pos++;
        }
    } else {
        id[HAtom(atype, ids[word.substr(from, to-from)])] = me;
    }
    return me;
}

double CherryForest::cherryness(const HAtom& a, const HAtom& b) {
    assert(a.atype() == b.atype());
    if (sizes[a.atype()] < 4) return 1.0;
    double positive = 0., total = 0.;
    for(CherryTree& ct : trees[a.atype()]) {
        positive += double(ct.has_cherry(a,b)) * ct.probability;
        total += ct.probability;
    }
    // if maximum probability is low, we cant trust cherryness much, 
    // therefore we add some bias to deal with it
    if (total < epsilon) total = 1;
    return min(1.0, positive / total + bias[a.atype()]);
}

void CherryForest::merge(const HAtom& a, const HAtom& b) {
    assert(a.atype() == b.atype());
    if (sizes[a.atype()] < 4) return;
    vector<CherryTree> new_trees;
    for(CherryTree& ct : trees[a.atype()]) {
        if (!ct.has_cherry(a,b)) continue;
        ct.merge(a,b);
        new_trees.push_back(move(ct));
    }
    trees[a.atype()] = move(new_trees);
    sizes[a.atype()]--;
}
    
void CherryForest::read_atom(History* history, int type, string filename) {
    trees[type].clear();
    sizes[type] = 0;
    for(auto la : history->leaf_atoms) for (auto atom : la.second) 
        sizes[type] += (atom.atype() == type);
    if (sizes[type] < 4) return;
    if (debuging)
        cout << "reading starts " << type << " " << filename << " " << sizes[type] << endl;
    ifstream file (filename, fstream::in);
    if (!file.is_open()) {
        cout << "ERROR: Cherry files not found" << endl;
        error_happened = 1;
        return;
    }

    string line, first, second, word;
    while(getline(file, line)) if (line == "   translate") break;
    map<string, int> ids;
    For(i, sizes[type]) {
        file >> first >> second;
        second.pop_back();
        ids[first] = HAtom::str_to_id(second);
    }
    //tree tree_1 [p = 0.780, P = 0.780] = [&W 0.780293] ((4,2),3,1);
    double probability = 0.0, best_probability = 0.0, cnt = 0.;
    while(file >> word) {
        if (word == "end;") break;
        assert(word == "tree");
        file >> word >> word >> word >> word >> word >> word >> word >> word >> word;
        file >> probability >> word >> word;
        word.pop_back();
        trees[type].push_back(CherryTree(probability, type, ids, word));
        best_probability = max(best_probability, probability);
        cnt++;
    }
    bias[type] = pow(1.-pow(best_probability,0.5),2.0);
    file.close();
}
    
CherryForest::CherryForest() {

}

CherryForest::CherryForest(CherryForest* original) {
    if (original == nullptr) return;
    this->trees = original->trees;
    this->sizes = original->sizes;
    this->bias = original->bias;
}
