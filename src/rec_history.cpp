#include <cassert>

#include "history.h"
#include "random.h"
#include "dynamics.h"
#include "score.h"

using namespace std;

void History::rec_parent(HEvent* event, Candidate* dont_use) {
    set<Candidate> cset = rec_candidates(event);
    vector<Candidate> cs(cset.begin(), cset.end());
    if (cs.size() == 0) return;
    vector<double> scores;
    vector<double> sum_scores = {0};

    if (neighbor_selection == neighbor_selection_enum::prioritize_used_events) {
        vector<bool> already_used;
        double used_scores_sum = 0;
        double unused_scores_sum = 0;
        for (auto c : cs) {
            double score = rec_score(c, event);
            scores.push_back(score);
            bool is_already_used = machine->was_duplication_used(c, event);
            already_used.push_back(is_already_used);
            if (is_already_used) {
                used_scores_sum += score;
            } else {
                unused_scores_sum += score;
            }
        }

        double current_ratio = (used_scores_sum) / (used_scores_sum + unused_scores_sum);
        double desired_ratio = prob_previously_used_event;
        if (current_ratio < desired_ratio) {
            double increase_ratio = desired_ratio / current_ratio;

            for (size_t i = 0; i < scores.size(); i++) {
                if (already_used[i]) {
                    scores[i] *= increase_ratio;
                }
            }
        }
    } else {
        for (auto c : cs) {
            if ((dont_use != nullptr) && (*dont_use == c)) {
                continue;
            }
            double score = rec_score(c, event);
            scores.push_back(score);
        }
    }

    for (auto score : scores) {
        sum_scores.push_back(sum_scores.back() + score);
    }

    double pick = random_double(0, sum_scores.back());
    For(i, cs.size()) if (pick < sum_scores[i+1]) {
        apply_candidate(cs[i], event);
        return;
    }
}

void History::apply_candidate(const Candidate& c, HEvent* event) {
    used_candidates.push_back(c);
    machine->add_used_duplication(c, event);
    rec_compute_parent(c, event);
    rec_merge_candidate(c, event);
}

set<Candidate> History::rec_candidates(HEvent* event) {
    set<Candidate> res;
    if (strategy == KNOW_HOW) {
        Dynamics d0 = Dynamics(this, event);
        d0.compute_graph(3.,0.9,0.4);
        Candidate c = d0.get_candidate(true);
        if (c.is_valid()) res.insert(c);
        return res;
    }
    if (strategy == NO_STRATEGY || strategy == CHERRY_NO ||
        strategy == CHERRY_TREE || strategy == CHERRY_LEN) {
        Dynamics d = Dynamics(this, event);
        d.compute_graph(2.0,0.5,0.05);
        int cnt = 0;
        while(true) {
            Candidate c = d.get_candidate();
            if (c.is_valid()) {
                res.insert(c);
                break;
            }
            assert(cnt++<1000);
        }
        return res;
    }

    Dynamics d = Dynamics(this, event);
    d.compute_graph(1.2,0.5,0.05);
    For(i, event->atoms.size()) {
        Candidate c = d.get_candidate();
        if (c.is_valid()) res.insert(c);
    }
    d.compute_graph(2.0,0.5,0.05);
    For(i, (event->atoms.size() * 3)) {
        Candidate c = d.get_candidate();
        if (c.is_valid()) res.insert(c);
    }
    for(auto c : res) {
        c.swap_dir();
        res.insert(c);
    }

    return res;
}

double History::rec_score(const Candidate& c, HEvent* event) {
    if (machine == nullptr) return 1.0;
    vdo values = all_scores(this, c, event);
    if (strategy == KNOW_HOW) {
        //machine->train_data(values, 1);
        return 1.0;
    }
    double res = machine->predict(values);
    return res;
}

void History::proc_learn() {
    HEvent* current = events.begin()->second;
    double now_time = current->event_time;
    current = new HEvent(current->species, gen_event_name(), "", current);
    events[current->name] = current;
    current->event_time = now_time -= 0.01;

    while(!current->is_final()) {
        set<Candidate> cs = rec_candidates(current);
        vector<vdo> good_values;
        vector<vdo> bad_values;
        for(auto c : cs) {
            vdo values = all_scores(this, c, current);
            HEvent *e = rec_see_event(c,current);
            if (is_original(e, strict_compare))
                good_values.push_back(values);
            else {
                if (strict_compare != SPECIAL_TRAINING || !is_original(e))
                    bad_values.push_back(values);
            }
            delete e;
            stats["candidates"] += 1;
        }
        For(i, min(good_values.size(), bad_values.size()))
            machine->train_data(good_values[i],1);
        For(i, min(good_values.size(), bad_values.size()))
            machine->train_data(bad_values[i],0);


        int old_strategy = strategy;
        int old_mode = cherry_mode;
        strategy = KNOW_HOW;
        cherry_mode = KNOW_HOW;
        rec_parent(current);
        strategy = old_strategy;
        cherry_mode = old_mode;

        if (current->parent == nullptr) return;
        bool correct = is_correct(true);
        while(current->parent != nullptr) {
            current = current->parent;
            events[current->name] = current;
            current->event_time = now_time -= 0.01;
        }
        if (!correct) return;
    }
    if (is_correct()) stats["full"] = 1;
}


void History::proc_reconstruct(int number) {
    // only linear reconstruction implemented
    assert(leaf_species.size()==1);
    // find last event
    HEvent* current = events.begin()->second;
    for(auto e : events)
        if (current->event_time > e.second->event_time)
            current = e.second;
    double now_time = current->event_time;
    if (current->type != "") {
        current = new HEvent(current->species, gen_event_name(), "", current);
        events[current->name] = current;
        current->event_time = now_time -= 0.01;
    }
    // reconstruct
    while(number!=0 && !current->is_final()) {
        if (number>0) number--;
        rec_parent(current);
        if (current->parent == nullptr) {
            assert(strategy == KNOW_HOW);
            return;
        }
        if (is_correct(true)) stats["max_events"] = events.size()-1;
        else if (number == EVAL_LAZY) return;
        while(current->parent != nullptr) {
            current = current->parent;
            events[current->name] = current;
            current->event_time = now_time -= 0.01;
        }
    }
}

bool History::real_reconstruct(HEvent* current) {
    double now_time;
    if (current == nullptr) { // start reconstruction from leafs
        current = events.begin()->second;
        now_time = current->event_time;
        current = new HEvent(current->species, gen_event_name(), "", current);
        events[current->name] = current;
        current->event_time = now_time -= 0.01;
    } else { // continue reconstruction from event
        now_time = current->event_time;
    }

    while(!current->is_final()) {
        rec_parent(current);
        if (current->parent == nullptr) {
            // TODO: This should never happen. You can always remove duplicated atoms, even one by one.
            cout << "Unsuccessful reconstruction" << endl;
            return false;
        }
        while(current->parent != nullptr) {
            current = current->parent;
            events[current->name] = current;
            current->event_time = now_time -= 0.01;
        }
    }

    machine->reset_used_duplications();
    return true;
}

void History::proc_test_candi(int strategy, string mark) {
    assert(leaf_species.size()==1);
    HEvent* current = leaf_events.begin()->second;
    set_strategy(strategy);

    int tot_cnt = 1000;
    int ok_cnt = 0;
    int oke_cnt = 0;
    int size_sum = 0;

    For(i, tot_cnt) {
        set<Candidate> cs = rec_candidates(current);
        assert(cs.size());
        int best = 0;
        for(auto c : cs) {
            HEvent *e = rec_see_event(c,current);
            best = max(best, is_original(e, strict_compare));
            delete e;
        }
        ok_cnt += bool(best);
        oke_cnt += (best==SAME_LAST);
        size_sum += cs.size();
    }

    stats["_size "+mark] = double(size_sum)/tot_cnt;
    stats["etest "+mark] = double(oke_cnt)/tot_cnt;
    stats["ctest "+mark] = double(ok_cnt)/tot_cnt;
}


void History::proc_test_score(int strategy, string mark) {
    assert(leaf_species.size()==1);
    HEvent* current = leaf_events.begin()->second;
    Machine* machine = nullptr;
    if (strategy == SCORE_CL) machine = new MachineOne();
    if (strategy == SCORE_BAC_NC) machine = new MachineBachelor();
    if (strategy == SCORE_BAC) machine = new MachineBachelor();
    if (strategy == SCORE_LR) machine = new MachineLinear();
    if (strategy == SCORE_LRS) machine = new MachineLinearStrict();
    if (machine != nullptr) machine->load();
    set_strategy(strategy, machine);

    int tot_cnt = 100;
    double sum_prob = 0.0;
    double sum_cnt = 0.0;
    int size_sum = 0;
    int is_max_sum = 0;

    For(i, tot_cnt) {
        double tot_score = 0.0;
        double ok_score = 0.0;
        double ok_cnt = 0;
        double max_score = 0.0;
        int is_max_good = 0;
        set<Candidate> cs = rec_candidates(current);
        assert(cs.size());

        for(auto c : cs) {
            double score = rec_score(c,current);
            HEvent *e = rec_see_event(c,current);
            if (score > max_score) is_max_good = 0;
            if(is_original(e, strict_compare)) {
                if (score > max_score) is_max_good = 1;
                ok_score += score;
                ok_cnt += 1;
            }
            max_score = max(max_score, score);
            //if (i == 0) {
            //    cout << mark << "score " << is_original(e) << " is " << score << endl;
            //}
            tot_score += score;
            delete e;
        }
        is_max_sum += is_max_good;
        sum_prob += ok_score / tot_score;
        sum_cnt += ok_cnt / double(cs.size());
        size_sum += cs.size();
    }
    if (machine != nullptr) delete machine;

    stats["_size "+mark] = double(size_sum)/tot_cnt;
    stats["cnt   "+mark] = double(sum_cnt)/tot_cnt;
    stats["score "+mark] = double(sum_prob)/tot_cnt;
    stats["good  "+mark] = double(is_max_sum)/tot_cnt;
}

History* History::rec_similar(int strategy, Machine* machine) {
    int last_used = random_int(0, used_candidates.size());

    History* h_new = new History(this);
    h_new->set_strategy(strategy, machine);

    HEvent* event = h_new->leaf_events["unicorn"];
    double now_time = event->event_time;
    event = new HEvent(event->species, h_new->gen_event_name(), "", event);
    h_new->events[event->name] = event;
    event->event_time = now_time -= 0.01;

    for (int i = 0; i < last_used; i++) {
        h_new->apply_candidate(used_candidates[i], event);
        if (event->parent == nullptr) {
            // TODO: This should never happen. You can always remove duplicated atoms, even one by one.
            cout << "Unsuccessful reconstruction" << endl;
            delete h_new;
            return nullptr;
        }
        while(event->parent != nullptr) {
            event = event->parent;
            h_new->events[event->name] = event;
            event->event_time = now_time -= 0.01;
        }
    }

    h_new->rec_parent(event, &used_candidates[last_used]);
    if (event->parent == nullptr) {
        // TODO: This should never happen. You can always remove duplicated atoms, even one by one.
        cout << "Unsuccessful reconstruction" << endl;
        delete h_new;
        return nullptr;
    }
    while(event->parent != nullptr) {
        event = event->parent;
        h_new->events[event->name] = event;
        event->event_time = now_time -= 0.01;
    }

    h_new->real_reconstruct(event);

    return h_new;
}
