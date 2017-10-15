#include"history.h"

void History::rec_parent(HEvent* event) {
    set<Candidate> cset = rec_candidates(event);
    vector<Candidate> cs(cset.begin(), cset.end());
    if (SIZE(cs) == 0) return;
    vector<double> sum_scores = {0};
    for(auto c : cs) {
        double score = rec_score(c, event);
        sum_scores.push_back(sum_scores.back() + score);
    }
    double pick = random_double(0, sum_scores.back());
    For(i, SIZE(cs)) if (pick < sum_scores[i+1]) {
        rec_compute_parent(cs[i], event);
        rec_merge_candidate(cs[i], event);
        return;
    }
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
    For(i, SIZE(event->atoms)) {
        Candidate c = d.get_candidate();
        if (c.is_valid()) res.insert(c);
    }
    d.compute_graph(2.0,0.5,0.05);
    For(i, (SIZE(event->atoms) * 3)) {
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
        For(i, min(SIZE(good_values), SIZE(bad_values))) 
            machine->train_data(good_values[i],1);
        For(i, min(SIZE(good_values), SIZE(bad_values)))
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
    assert(SIZE(leaf_species)==1);
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
        if (is_correct(true)) stats["max_events"] = SIZE(events)-1;
        else if (number == EVAL_LAZY) return;
        while(current->parent != nullptr) {
            current = current->parent;
            events[current->name] = current;
            current->event_time = now_time -= 0.01;
        }
    }
}

void History::real_reconstruct() {
    HEvent* current = events.begin()->second;
    double now_time = current->event_time;
    current = new HEvent(current->species, gen_event_name(), "", current);
    events[current->name] = current;
    current->event_time = now_time -= 0.01;
    
    while(!current->is_final()) {
        rec_parent(current);
        if (current->parent == nullptr) {
            cout << "Unsuccessfull reconstruction" << endl;
            exit(1);
        }
        while(current->parent != nullptr) {
            current = current->parent;
            events[current->name] = current;
            current->event_time = now_time -= 0.01;
        }
    }
}

void History::proc_test_candi(int strategy, string mark) {
    assert(SIZE(leaf_species)==1);
    HEvent* current = leaf_events.begin()->second;
    set_strategy(strategy);

    int tot_cnt = 1000;
    int ok_cnt = 0;
    int oke_cnt = 0;
    int size_sum = 0;

    For(i, tot_cnt) {
        set<Candidate> cs = rec_candidates(current);
        assert(SIZE(cs));
        int best = 0;
        for(auto c : cs) {
            HEvent *e = rec_see_event(c,current);
            best = max(best, is_original(e, strict_compare));
            delete e;
        }
        ok_cnt += bool(best);
        oke_cnt += (best==SAME_LAST);
        size_sum += SIZE(cs);
    }
    
    stats["_size "+mark] = double(size_sum)/tot_cnt;
    stats["etest "+mark] = double(oke_cnt)/tot_cnt;
    stats["ctest "+mark] = double(ok_cnt)/tot_cnt;
}


void History::proc_test_score(int strategy, string mark) {
    assert(SIZE(leaf_species)==1);
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
        assert(SIZE(cs));

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
        sum_cnt += ok_cnt / double(SIZE(cs));
        size_sum += SIZE(cs);
    }
    if (machine != nullptr) delete machine;
    
    stats["_size "+mark] = double(size_sum)/tot_cnt;
    stats["cnt   "+mark] = double(sum_cnt)/tot_cnt;
    stats["score "+mark] = double(sum_prob)/tot_cnt;
    stats["good  "+mark] = double(is_max_sum)/tot_cnt;
}

