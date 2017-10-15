#include"gatom.h"

string invert_dna(const string& str, bool really = true) {
    if (!really) return str;
    string res = string(SIZE(str), 0);
    For(i, SIZE(str)) res[i] = base_inv[int(str[SIZE(str)-1-i])];
    return res;
}

// remove all '-'
string compact_dna(const string& str, bool really = true) {
    if (!really) return str;
    string res = "";
    for(auto c : str) if (c!='-') res.push_back(c);
    return res;
}

int GAtomType::type_cnt = 0;

void GAtomType::split(int position) {
    this->root->split(position, nullptr, nullptr);
}
void GAtomType::invert() {
    this->root->hard_invert();
}
GAtomType::GAtomType(GAtom* root) {
    this->root = root;
    this->id = ++type_cnt;
}

bool GAtom::change_parent(GAtom* new_parent) {
    if (new_parent != nullptr) this->type = new_parent->type;
    if (parent == new_parent) return false;
    if (parent != nullptr) { 
        auto it = remove(parent->children.begin(), parent->children.end(), this);
        parent->children.resize(distance(parent->children.begin(), it));
    }
    if (new_parent != nullptr)
        new_parent->children.push_back(this);
    parent = new_parent;
    return true;
}

void GAtom::invert() {
    inverted = !inverted; 
}
void GAtom::hard_invert() {
    inverted = !inverted;
    dna = invert_dna(dna);
    for(GAtom* child : children)
        child->hard_invert();
}

void GAtom::mutate(double time) {
    for(auto& c : dna)
        if (c != '-')
            c = Model::instance()->get_mutated_base(c, time);
    For(i, SIZE(dna))
        if (Model::instance()->get_indel_happened(time)) {
            do {
                dna[i] = '-';
                i++;
            } while(i < SIZE(dna) && rand()%2 == 0);
        }
    
}
void GAtom::split(int position, GAtom* first_parent, GAtom* second_parent) {
    if (inverted) {
        GAtom* first = new GAtom(first_parent, dna.substr(0, position));
        first->inverted = inverted;
        dna = dna.substr(position);
        this->change_parent(second_parent); 
        first->next = this->next;
        this->next = first;
        vector<GAtom*> children_copy = children;
        for(GAtom* child : children_copy)
            child->split(position, first, this);
    } else {
        GAtom* second = new GAtom(second_parent, dna.substr(position));
        second->inverted = inverted;
        dna.resize(position);
        this->change_parent(first_parent);
        second->next = this->next;
        this->next = second;
        vector<GAtom*> children_copy = children;
        for(GAtom* child : children_copy)
            child->split(position, this, second);
    } 
}
    
GAtom::GAtom(int length) {
    this->inverted = 0;
    this->next = nullptr;
    this->parent = nullptr;
    this->type = new GAtomType(this);
    dna = "";
    For(i, length) dna += bases[rand()%BASES];
}

GAtom::GAtom(GAtom* parent, const string &dna) {
    if ((this->parent = parent) == nullptr) {
        this->inverted = 0;
        this->type = new GAtomType(this);
        this->next = nullptr;
    } else {
        this->type = parent->type;
        this->inverted = parent->inverted;
        this->next = parent->next;
        parent->children.push_back(this);
    }
    this->dna = dna;
}

void GAtom::write_dna(ostream& os, const string& sep, bool compact) {
    os << invert_dna(compact_dna(dna, compact), inverted) << sep;
}
void GAtom::write_type(ostream& os, const string& sep) {
    os << get_id() << sep;
}
ostream& operator<<(ostream& os, GAtom& atom) {
    return os << atom.get_id() << "(" << atom.length() << ")";
}
