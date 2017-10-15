#ifndef SCORE_H
#define SCORE_H

#include"constants.h"

#define SC_SEQ_LEN 0
#define SC_SEQ_NUM_TYPES 1
//#define SC_SEQ_DNA_LEN -1
// estimated events to end
#define SC_POST_BPC 2

#define SC_EV_LEN 3
#define SC_EV_DIST 4
#define SC_DEL_NUM 5
#define SC_DEL_LEN 6
#define SC_EV_SIDES 7
#define SC_EV_PREV_SP 8
#define SC_EV_POST_SP 9
#define SC_EV_PREV_BP 10
#define SC_EV_POST_BP 11
#define SC_EV_DSIG 12

#define SC_AVG_CHERRY 13
#define SC_PROD_CHERRY 14
#define SC_LEN_CHERRY 15

vector<double> all_scores(History* h, const Candidate& c, HEvent* e);

class ScoringData {
public:
    History* history;
    HEvent *first_event, *last_event;
    Candidate c;
    vector<pii> deletions;
    vector<pair<HAtom,HAtom>> atom_friends;
    set<pii> prev_bp, post_bp;
    set<pii> prev_sp, post_sp;
    int prev_bpc, post_bpc, num_types;
    //map<int,int> prev_leftnum, prev_rightnum, post_lefnum, post_rightnum;
    
    ScoringData(History* h, const Candidate& c, HEvent* e);
    ~ScoringData();
};


#endif
