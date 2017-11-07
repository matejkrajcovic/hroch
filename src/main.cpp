#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <cassert>

#include "ghistory.h"
#include "constants.h"

using namespace std;

void generate_history(double time, string id) {
    GHistory* history = new GHistory();
    history->generate_random(time, 100000);
    history->save_to_files(datapath + "generated", id);
    delete history;
}

void generate_test() {
    for(int i = 100; i < 200; ++i) generate_history(0.02, "F2"+to_string(i));
    for(int i = 100; i < 200; ++i) generate_history(0.03, "F3"+to_string(i));
    for(int i = 100; i < 200; ++i) generate_history(0.04, "F4"+to_string(i));
}

void generate_all() {
    // if random seed is same, we have to generate histories sequential
    generate_test();
    for(int i = 10000; i < 15000; ++i) generate_history(0.04, "L4"+to_string(i));
}

void write_all_stats(map<string, vector<double>>& stats, ostream& os) {
    os << "all stats" << endl;
    for(auto sp : stats) {
        vector<double> stat = sp.second;
        sort(stat.begin(), stat.end());
        double avg = 0;
        for(double d : stat) avg += d;
        avg /= stat.size();
        double p0 = stat[0];
        double p10 = stat[int(stat.size()*0.1)];
        double p50 = stat[int(stat.size()*0.5)];
        double p90 = stat[int(stat.size()*0.9)];

        os << sp.first << "  AGV: " << avg
           << "   0% " <<  p0 << "  10% " << p10
           << "  50% " << p50 << "  90% " << p90 << endl;
    }
}

void test_candi() {
    map<string, vector<double>> stat_val;
    string prefix = TEST_CASE;
    ofstream file("stats/candi-stats"+prefix);

    for(int hid = LOWER_RANGE; hid<UPPER_RANGE; ++hid) {
        History* h0 = new History(datapath + "generated", prefix+to_string(hid));
        History* h1 = new History(h0);
        h1->stats["_name"] = hid;

        h1->proc_test_candi(CHERRY_NO,"1nc");
        h1->proc_test_candi(CHERRY_TREE,"1ct");
        h1->proc_test_candi(CHERRY_LEN,"1cl");
        h1->proc_test_candi(SCORE_BAC_NC,"nc ");
        h1->proc_test_candi(SCORE_BAC,"def");
        h1->proc_test_candi(SCORE_CL,"cl ");

        h1->write_stats(file);
        file << endl;
        for(auto sp : h1->stats) stat_val[sp.first].push_back(sp.second);
        delete h0;
        delete h1;
    }
    write_all_stats(stat_val,cout);
    write_all_stats(stat_val,file);
    file.close();
}

void test_score() {
    map<string, vector<double>> stat_val;
    string prefix = TEST_CASE;
    string strictsuffix = "";
    if (strict_compare) strictsuffix = "strict";
    ofstream file("stats/score-stats"+prefix+strictsuffix);

    for(int hid = LOWER_RANGE; hid<UPPER_RANGE; ++hid) {
        History* h0 = new History(datapath + "generated", prefix+to_string(hid));
        History* h1 = new History(h0);
        h1->stats["_name"] = hid;

        h1->proc_test_score(CHERRY_TREE,"1ct");
        h1->proc_test_score(SCORE_BAC_NC,"nc ");
        h1->proc_test_score(SCORE_BAC,"def");
        h1->proc_test_score(SCORE_LR,"lr ");
        h1->proc_test_score(SCORE_LRS,"lrs");

        h1->write_stats(file);
        file << endl;
        for(auto sp : h1->stats) stat_val[sp.first].push_back(sp.second);
        delete h0;
        delete h1;
    }
    write_all_stats(stat_val,cout);
    write_all_stats(stat_val,file);
    file.close();
}

void train_history(Machine* machine, string name) {
    error_happened = 0;
    History* h0 = new History(datapath + "generated", name);
    if (error_happened) {
        cout << "Could not train with hitory " << name << endl;
        return;
    }

    History* hsp = new History(h0);
    hsp->set_strategy(SCORE_LR,machine);
    hsp->proc_learn();
    //hsp->save("g-train-"+name);

    delete hsp;
    delete h0;
}

void train() {
    Machine* machine = new MachineLinear();
    for(int i = 10000; i < 15000; ++i) train_history(machine,"L4"+to_string(i));
    machine->save();
    delete machine;
    cout << "training finished" << endl;
}

void reconstruct_one(History* h0, string hid, int strategy) {
    cout << "Reconstructing " << hid << " with strategy " << strategy << endl;

    Machine* machine = nullptr;
    if (strategy == SCORE_CL) machine = new MachineOne();
    if (strategy == SCORE_BAC_NC) machine = new MachineBachelor();
    if (strategy == SCORE_BAC) machine = new MachineBachelor();
    if (strategy == SCORE_LR) machine = new MachineLinear();
    if (strategy == SCORE_LRS) machine = new MachineLinearStrict();
    if (machine != nullptr) machine->load();

    int attempts = 1000;
    int max_events = 0;
    double avg_num_events = 0.0;
    vector<int> dist_num_events(1000,0);
    vector<int> cnt_events(h0->events.size()-1,0);

    For(i, attempts) {
        if (i>0 && i%2000==0) cout << "history " << hid << "  attempt " << i << endl;
        History* h1 = new History(h0);
        h1->set_strategy(strategy, machine);
        h1->proc_reconstruct(stats?EVAL_INF:EVAL_LAZY);
        max_events = max(max_events, int(h1->stats["max_events"]));
        cnt_events[int(h1->stats["max_events"])]++;
        if (stats) {
            avg_num_events += h1->events.size();
            dist_num_events[h1->events.size()]++;
        }
        delete h1;
    }
    avg_num_events /= attempts;
    while(dist_num_events.size() && dist_num_events.back() == 0)
        dist_num_events.pop_back();
    dist_num_events[h0->events.size()-1]--;

    cout << "  " << h0->events.size() << " events,  "
         << h0->leaf_atoms.begin()->second.size() << " atoms" << endl;
    cout << "  max_events " << max_events << "    " << cnt_events << endl;
    if (stats) {
        cout << "  length of events: avg=" <<avg_num_events << "  ::  ";
        cout << dist_num_events << endl;
    }
    cout << endl;
    if (machine != nullptr) delete machine;
}

void reconstruct_many(string hid) {
    History* h0 = new History(datapath + "generated", hid);
    if (debuging) {
        for(auto ev : h0->events) ev.second->test_stats(h0, cout);
        for(auto ev : h0->events) ev.second->write_detailed(cout);
    }
    History* hsp = new History(h0);
    hsp->set_strategy(KNOW_HOW);
    hsp->proc_reconstruct();
    if (!hsp->is_correct()) {
        cout << hid << " cannot be reconstructed" << endl;
        //delete hsp;
        //return;
    }
    delete hsp;
    //reconstruct_one(h0, hid, CHERRY_NO);
    reconstruct_one(h0, hid, CHERRY_TREE);
    //reconstruct_one(h0, hid, CHERRY_LEN);
    reconstruct_one(h0, hid, SCORE_BAC_NC);
    reconstruct_one(h0, hid, SCORE_BAC);
    reconstruct_one(h0, hid, SCORE_LR);
    reconstruct_one(h0, hid, SCORE_LRS);
    delete h0;
}

void reconstruct(string atoms_file, string trees_dir, int count, int strategy) {
    string outputfile_name = "hroch_" + atoms_file;
    for(auto& c : outputfile_name) if (c=='/') c = '-';
    outputfile_name = "outputs/"+outputfile_name+".histories";
    ofstream ofile(outputfile_name, fstream::out);

    Machine* machine = nullptr;
    if (strategy == SCORE_CL) machine = new MachineOne();
    if (strategy == SCORE_BAC_NC) machine = new MachineBachelor();
    if (strategy == SCORE_BAC) machine = new MachineBachelor();
    if (strategy == SCORE_LR) machine = new MachineLinear();
    if (machine != nullptr) machine->load();

    int shortest = 123456789;
    vector<History*> hists;
    vector<int> lengths(1000, 0);

    For(i, count) {
        History* h = new History(atoms_file, trees_dir, strategy);
        h->set_strategy(strategy, machine);
        h->real_reconstruct();
        int num_events = h->events.size();
        lengths[num_events]++;
        if (num_events < shortest) {
            shortest = num_events;
            for(auto hp : hists) delete hp;
            hists.clear();
        }
        if (num_events == shortest) {
            bool is_there = false;
            for(auto hist : hists) if (h->same_as(hist)) is_there = true;
            if (!is_there) hists.push_back(h);
        }
        h->write_events(ofile);
        ofile << endl;
    }
    while(lengths.back() == 0) lengths.pop_back();

    for(auto h : hists) {
        h->write_events(cout);
        cout << endl;
    }
    cout << shortest << endl;
    cout << lengths << endl;
    ofile << shortest << endl;
    ofile << lengths << endl;
    if (machine != nullptr) delete machine;
    ofile.close();
}





void help(string name) {
    cout << "Usage: "<< name << "--solve ATOMFILE TREESDIR COUNT [STRATEGY]" << endl;
    cout << "       Reconstructs COUNT histories of sequence given by ATOMFILE." << endl;
    cout << "       Outputs some statistics, and set of shortest histories." << endl;
    cout << "       Use one of te following strategies:" << endl;
    cout << "          1: no cherryness, no scoring" << endl;
    cout << "          2: use cherryness, but no scoring" << endl;
    cout << "          5: simple scoring, no cherryness" << endl;
    cout << "          6: simple scoring, with cherryness" << endl;
    cout << "          7: advanced scoring, with cherryness" << endl;
    cout << "          8: [default] advanced strict scoring, with cherryness" << endl;

    cout << endl << "other (undocumented) options:" << endl;
    cout << "  --gen-test:      generate test data" << endl;
    cout << "  --gen-all:       generate test and train data" << endl;
    cout << "  --train:         produce lr-train file" << endl;
    cout << "  --test-c:        make statistics of candidate proposal algorithms" << endl;
    cout << "  --test-s:        make statistics of scoring algorithms" << endl;
    cout << "  --rec:           make statistics of history reconstructin algorithms" << endl;
}


int main(int argc, char **argv) {
    setup_constants();
    auto args = parse_arguments(argc, argv);
    cout << "HROCH: Heuristic reconstrucion of cluster histories." << endl;
    if (args.count("--help")) {
        help(argv[0]);
        return 0;
    }
    if (args.count("--gen-test")) generate_test();
    if (args.count("--gen-all")) generate_all();
    if (args.count("--train")) train();
    if (args.count("--test-c")) test_candi();
    if (args.count("--test-s")) test_score();
    if (args.count("--rec")) {
        for(int i = LOWER_RANGE; i<UPPER_RANGE; ++i)
            reconstruct_many(TEST_CASE+to_string(i));
    }
    if (args.count("--solve")) {
        assert(argc >= 5);
        int strategy = SCORE_LR;
        if (argc > 5) strategy = stoi(string(argv[5]));
        int count = stoi(string(argv[4]));
        string atoms_file = argv[2];
        string trees_dir = argv[3];
        reconstruct(atoms_file, trees_dir, count, strategy);
    }
    cout << "Done." << endl;
}
