#ifndef CANDIDATE_H
#define CANDIDATE_H

#include <algorithm>
#include <vector>
#include <ostream>

class Candidate {
public:
    int b1,e1,b2,e2;
    inline bool is_inv() const { return e2 < b2; }
    inline bool is1_left() const { return e1 < e2; }
    inline bool is1_right() const { return e2 < e1; }
    inline bool is_valid() const { return e1 < std::min(b2+1, e2) || std::max(b2-1, e2) < b1+1; }
    std::vector<int> directions;
    Candidate(int b1, int e1, int b2, int e2);
    Candidate(int b1, int e1, int b2, int e2, const std::vector<int>& directions);

    void swap_dir();

    friend bool operator==(const Candidate& c1, const Candidate& c2) {
        return (c1.b1 == c2.b1) &&
               (c1.e1 == c2.e1) &&
               (c1.b2 == c2.b2) &&
               (c1.e2 == c2.e2) &&
               (c1.directions == c2.directions);
    }
    friend bool operator<(const Candidate& c1, const Candidate& c2) {
        if (c1.b1 != c2.b1) return c1.b1 < c2.b1;
        if (c1.e1 != c2.e1) return c1.e1 < c2.e1;
        if (c1.b2 != c2.b2) return c1.b2 < c2.b2;
        if (c1.e2 != c2.e2) return c1.e2 < c2.e2;
        return c1.directions < c2.directions;
    }
};

std::ostream& operator<<(std::ostream& os, const Candidate& c);

#endif
