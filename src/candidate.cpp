#include "candidate.h"
#include "utils.h"
using namespace std;

Candidate::Candidate(int b1, int e1, int b2, int e2) {
    this->b1 = b1;
    this->e1 = e1;
    this->b2 = b2;
    this->e2 = e2;
}

Candidate::Candidate(int b1, int e1, int b2, int e2, const vector<int>& directions) :
    Candidate::Candidate(b1,e1,b2,e2) {
    this->directions = directions;
}

void Candidate::swap_dir() {
    if (is_inv()) {
        swap(b1,e2);
        swap(b2,e1);
        b1--;
        e1--;
        b2++;
        e2++;
        reverse(directions.begin(), directions.end());
    } else {
        swap(b1,b2);
        swap(e1,e2);
    }
    for(auto& d : directions) {
        if (d == 1) d = 2;
        else if (d == 2) d = 1;
    }
}

ostream& operator<<(ostream& os, const Candidate& c) {
    return os << "(" << c.b1 << " " << c.e1 << " -> "
       << c.b2 << " " << c.e2 << " " << c.directions << ")";
}
