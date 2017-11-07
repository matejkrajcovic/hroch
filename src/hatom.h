// Representation of atom used in reconstruction

#ifndef HATOM_H
#define HATOM_H

#include <map>
#include <vector>
#include <string>
#include <ostream>

class HAtom {
    static std::map<std::string, int> si_map;
    static std::map<int, std::string> is_map;

    std::vector<int> ids;
public:
    static std::string id_to_str(const int& id);
    static int str_to_id(const std::string& str);
    static void clear_strid_mapping();

    void add_id(int id);
    void add_ids(const std::vector<int>& what);
    std::vector<int> get_ids() const;
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
    friend std::ostream& operator<<(std::ostream& os, const HAtom& a) {
        return os << a.type;
    }
};

#endif
