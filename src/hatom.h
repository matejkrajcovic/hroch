// Representation of atom used in reconstruction

#ifndef HATOM_H
#define HATOM_H

#include"constants.h"

class HAtom {
    static map<string, int> si_map;
    static map<int, string> is_map;
    
    vector<int> ids;
public:
    static string id_to_str(const int& id);
    static int str_to_id(const string& str);
    static void clear_strid_mapping();

    void add_id(int id);
    void add_ids(const vector<int>& what);
    vector<int> get_ids() const;
    int type;
    int atype() const { return abs(type); }

    HAtom(int type);
    HAtom(int type, int id);

    friend bool operator==(const HAtom& h1, const HAtom& h2) {
        return h1.atype() == h2.atype() && h1.ids == h2.ids;
    }
    friend bool operator<(const HAtom& h1, const HAtom& h2) {
        return (h1.atype() != h2.atype())?(h1.atype() < h2.atype()):(h1.ids < h2.ids);
    }
    friend ostream& operator<<(ostream& os, const HAtom& a) {
        return os << a.type;
    }
};


#endif
