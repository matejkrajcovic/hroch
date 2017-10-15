#include"hatom.h"

map<string, int> HAtom::si_map;
map<int, string> HAtom::is_map;

string HAtom::id_to_str(const int& id) {
    assert(is_map.count(id));
    return is_map[id];
}

int HAtom::str_to_id(const string& str) {
    if (si_map.count(str) == 0) {
        si_map[str] = SIZE(si_map);
        is_map[si_map[str]] = str;
    }
    return si_map[str];
}

void HAtom::clear_strid_mapping() {
    is_map.clear();
    si_map.clear();
}

HAtom::HAtom(int type) {
    this->type = type;
}
HAtom::HAtom(int type, int id) : 
    HAtom::HAtom(type) {
    this->add_id(id);
}

void HAtom::add_id(int id) {
    For(i, SIZE(ids)) {
        if (ids[i] == id) return;
        if (ids[i] > id) swap(ids[i], id);
    }
    ids.push_back(id);
}
void HAtom::add_ids(const vector<int>& what) {
    ids.insert(ids.end(), what.begin(), what.end());
    sort(ids.begin(), ids.end());
    auto it = unique(ids.begin(), ids.end());
    ids.resize(distance(ids.begin(), it));
}

vector<int> HAtom::get_ids() const {
    return ids;
}
