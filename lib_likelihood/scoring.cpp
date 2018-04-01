#include"event.h"
#include"atomtree.h"

namespace likelihood {
double best_time(const vector<pair<Node*, Node*> >& nodes){
    double total = 0.0, same = 0.0;
    for(auto np : nodes) For(i, SIZE(np.first->allprobs[0])){
        double len1 = 0.0, len2 = 0.0;
        For(j, 4) {
            len1 += (np.first->allprobs[j][i])*(np.first->allprobs[j][i]);
            len2 += (np.second->allprobs[j][i])*(np.second->allprobs[j][i]);
            same += (np.first->allprobs[j][i])*(np.second->allprobs[j][i]);
        }
        total += sqrt(len1*len2);
    }
    assert(total > 0.01 && same > total*0.3);
    //cerr << same << " " << total << " " << (4*same/(3*total)) << " " << -3./4. * log((4.*same)/(3.*total) -1./3.) << endl;
    return -3./4. * log((4.*same)/(3.*total) -1./3.);

}

double Event::path_score(const vtri& path){
    vector<piii> deletions;
    vector<int> data;
    Event* e = new Event();
    compute_prevgenes(path, e, deletions, data);

    // Scoring length

    double sc_length = 1.0;
    if (history->scoring_length){
        for(auto x : path) sc_length += double(x.a%3==0);
    }
    sc_length = log(sc_length)*history->scoring_length;

    // Scoring oldhistory

    //int &adi1 = data[1], &adi2 = data[2], &adi3 = data[3]; //, &adi4 = data[4], &adi5 = data[5];
    double sc_oldhistory = 0;
    if (history->scoring_oldhistory){
        //vi change = vi(e->genes.begin()+adi1+adi2, e->genes.begin()+adi1+adi2+adi3);

        auto it = history->old_events.begin();
        //if ((it = history->old_events.find(change)) != history->old_events.end())
        //    sc_oldhistory += it->second;
        if ((it = history->old_seqences.find(e->genes)) != history->old_seqences.end())
            sc_oldhistory += it->second;
        sc_oldhistory*=log(1.0+sc_oldhistory)*history->scoring_oldhistory;
    }

    // Scoring boundary

    double sc_boundary = 0.0;
    if (history->scoring_boundary){
        double type_bonus = 0.0;
        if (path.front().a >=3){
            set<int> types;
            for(auto x : genes) types.insert(abs(x));
            int typediff = SIZE(set<int>(genes.begin(), genes.end())) - SIZE(set<int>(e->genes.begin(), e->genes.end()));
            int typeleft = SIZE(genes) - SIZE(types);
            type_bonus += log(1.0+typediff/(1.0+typeleft)) - 0.5;
        }
        set<pii> typepairs;

        For(i, SIZE(e->genes)+1) {
            int x = (i==0)?0:e->genes[i-1];
            int y = (i==SIZE(e->genes))?0:e->genes[i];
            typepairs.insert(pii(x,y));
            typepairs.insert(pii(-y,-x));
        }
        For(i, SIZE(genes)+1) {
            int x = (i==0)?0:genes[i-1];
            int y = (i==SIZE(genes))?0:genes[i];
            if (typepairs.count(pii(x,y))) continue;
            type_bonus+=1;
        }

        sc_boundary = history->scoring_boundary*type_bonus;
    }

    // Scoring similarity

    double gap_time = 0.01;
    double std_deviation = 0.0;
    double sc_similar = 0;
    if (history->scoring_similar){
        int totallength = 0;
        vector<pair<Node*, Node*> > nodepairs;
        for(auto x : path) if (x.a%3==0) {
            assert(SIZE(atom_nodes[x.b]->allprobs[0])==SIZE(atom_nodes[x.c]->allprobs[0]));
            nodepairs.push_back(make_pair(atom_nodes[x.b], atom_nodes[x.c]));
            totallength += SIZE(atom_nodes[x.b]->allprobs[0]);
        }
        gap_time = best_time(nodepairs);
        //cerr << "DPTNG " << SIZE(deletions) << " " << SIZE(path) << " ";
        //cerr << totallength << " " << SIZE(nodepairs) << " " << gap_time << endl;
        for(auto x : path) if(x.a%3==0){
            nodepairs.clear();
            nodepairs.push_back(make_pair(atom_nodes[x.b], atom_nodes[x.c]));
            double part_gap_time = best_time(nodepairs);
            int length = SIZE(atom_nodes[x.b]->allprobs[0]);
            std_deviation += length * (part_gap_time-gap_time)*(part_gap_time-gap_time);
        }
        std_deviation = sqrt(std_deviation/totallength);
        sc_similar = log(1.0+double(totallength)/(gap_time+0.01))*
                     log(1.0+double(totallength)/(std_deviation+0.01))*0.1*history->scoring_similar;

    }

    // Scoring cherryness
    double sc_cherryness = 0.0;

    // Conclusion
    double score = sc_length + sc_oldhistory + sc_boundary + sc_similar + sc_cherryness - 8.0;
    //fprintf(stderr, "aaa %lf %lf %lf %lf %lf %lf\n",score, sc_length,
    //        sc_oldhistory, sc_boundary, sc_similar, sc_cherryness);

    return exp(score/double(history->scoring_divisor));
}
}
